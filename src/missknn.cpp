// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <progress.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

// RcppProgress's master-thread detection (InterruptableProgressMonitor::is_master())
// is only correct under OpenMP: without _OPENMP defined it always returns true, so
// every RcppParallel worker thread would believe itself to be the master and race on
// the shared progress state. This package parallelizes via RcppParallel (TBB/std::thread),
// not OpenMP, so Progress::increment() must never be called from inside a Worker's
// operator(). When display_progress is requested, the same Worker functor is instead
// invoked serially, one chunk at a time, from the calling (single) thread, which is
// the documented-safe way to drive a Progress bar.
template <typename WorkerT>
static void missknn_run_with_optional_progress(WorkerT& worker, int n, bool display_progress) {
  if (!display_progress) {
    RcppParallel::parallelFor(0, n, worker);
    return;
  }
  Progress p(static_cast<unsigned long>(n), true);
  for (int i = 0; i < n; ++i) {
    if (Progress::check_abort()) {
      Rcpp::stop("missknn: interrupted by user.");
    }
    worker(static_cast<std::size_t>(i), static_cast<std::size_t>(i + 1));
    p.increment();
  }
}
using namespace RcppParallel;

// Minimum relative holdout-error improvement a local (KNN/regression)
// candidate must show over the global mean/mode baseline before tuning
// abandons the baseline. Guards against chance sampling noise in the
// holdout evaluation from sending a no-signal column down the slower,
// noisier local path for no real accuracy benefit.
static const double kTuneMargin = 0.08;

// Regression's chosen bandwidth is picked as the minimum holdout error over
// several candidates (see cpp_missknn_tune_numeric); taking a minimum over
// multiple noisy candidates is itself biased low relative to any single
// candidate's true error (the same effect that makes "best of several dice
// rolls" beat a fair average even when the die is unweighted), so accepting
// regression needs a stricter margin than a single-candidate comparison
// (like the plain weighted-mean k search) does.
static const double kRegressionMargin = 0.20;

// All hot-loop helpers below take RcppParallel's RMatrix/RVector wrappers
// (thread-safe pointer access, no R API calls) rather than Rcpp types
// directly, so the same code path can run both serially and inside a
// RcppParallel::parallelFor worker without a separate implementation.
static inline double missknn_distance(
    int receiver,
    int donor,
    int target,
    const RMatrix<double>& numeric_scaled,
    const RMatrix<int>& categorical_codes,
    const std::vector<int>& col_types,
    const std::vector<int>& numeric_pos,
    const std::vector<int>& categorical_pos,
    const std::vector<double>& weights) {
  double sum_w = 0.0;
  double sum_d = 0.0;
  const int p = static_cast<int>(col_types.size());

  for (int col = 0; col < p; ++col) {
    if (col == target) {
      continue;
    }

    const double w = weights[col];
    if (!(w > 0.0)) {
      continue;
    }

    if (col_types[col] == 1) {
      const int pos = numeric_pos[col] - 1;
      if (pos < 0) {
        continue;
      }
      const double r = numeric_scaled(receiver, pos);
      const double d = numeric_scaled(donor, pos);
      if (!R_IsNA(r) && !R_IsNA(d)) {
        const double diff = d - r;
        sum_d += w * diff * diff;
        sum_w += w;
      }
    } else {
      const int pos = categorical_pos[col] - 1;
      if (pos < 0) {
        continue;
      }
      const int r = categorical_codes(receiver, pos);
      const int d = categorical_codes(donor, pos);
      if (r != NA_INTEGER && d != NA_INTEGER) {
        sum_d += w * (r != d ? 1.0 : 0.0);
        sum_w += w;
      }
    }
  }

  if (sum_w <= 0.0) {
    return R_PosInf;
  }
  return std::sqrt(sum_d / sum_w);
}

// Partial top-k selection (average O(n)) instead of a full sort.
static void missknn_topk(
    const std::vector<double>& distances,
    const std::vector<int>& donor_rows,
    int k,
    std::vector<int>& out_idx) {
  const int n = static_cast<int>(distances.size());
  std::vector<int> finite_idx;
  finite_idx.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (R_finite(distances[i])) {
      finite_idx.push_back(i);
    }
  }

  const int kk = std::min(k, static_cast<int>(finite_idx.size()));
  out_idx.clear();
  if (kk <= 0) {
    return;
  }

  auto cmp = [&](int a, int b) {
    if (distances[a] == distances[b]) {
      return donor_rows[a] < donor_rows[b];
    }
    return distances[a] < distances[b];
  };

  std::nth_element(finite_idx.begin(), finite_idx.begin() + (kk - 1), finite_idx.end(), cmp);
  std::sort(finite_idx.begin(), finite_idx.begin() + kk, cmp);
  out_idx.assign(finite_idx.begin(), finite_idx.begin() + kk);
}

static double missknn_numeric_fallback(const RMatrix<double>& numeric_raw, const std::vector<int>& donor_rows, int target_pos) {
  double sum = 0.0;
  int n = 0;
  for (size_t i = 0; i < donor_rows.size(); ++i) {
    const double v = numeric_raw(donor_rows[i] - 1, target_pos);
    if (!R_IsNA(v)) {
      sum += v;
      ++n;
    }
  }
  return n > 0 ? sum / static_cast<double>(n) : NA_REAL;
}

static int missknn_categorical_fallback(const RMatrix<int>& categorical_codes, const std::vector<int>& donor_rows, int target_pos) {
  std::map<int, int> counts;
  for (size_t i = 0; i < donor_rows.size(); ++i) {
    const int v = categorical_codes(donor_rows[i] - 1, target_pos);
    if (v != NA_INTEGER) {
      counts[v] += 1;
    }
  }
  if (counts.empty()) {
    return NA_INTEGER;
  }
  int best_code = NA_INTEGER;
  int best_count = -1;
  for (const auto& kv : counts) {
    if (kv.second > best_count || (kv.second == best_count && kv.first < best_code)) {
      best_code = kv.first;
      best_count = kv.second;
    }
  }
  return best_code;
}

// Solve (X'WX + lambda I) b = X'Wy via Gauss-Jordan elimination.
// A is (q x q) row-major, rhs is length q. Solved in place; returns false if singular.
static bool solve_linear_system(std::vector<double>& A, std::vector<double>& rhs, int q) {
  for (int col = 0; col < q; ++col) {
    int pivot = col;
    double best = std::fabs(A[col * q + col]);
    for (int r = col + 1; r < q; ++r) {
      const double v = std::fabs(A[r * q + col]);
      if (v > best) {
        best = v;
        pivot = r;
      }
    }
    if (best < 1e-10) {
      return false;
    }
    if (pivot != col) {
      for (int c = 0; c < q; ++c) {
        std::swap(A[col * q + c], A[pivot * q + c]);
      }
      std::swap(rhs[col], rhs[pivot]);
    }
    const double diag = A[col * q + col];
    for (int c = 0; c < q; ++c) {
      A[col * q + c] /= diag;
    }
    rhs[col] /= diag;
    for (int r = 0; r < q; ++r) {
      if (r == col) continue;
      const double factor = A[r * q + col];
      if (factor == 0.0) continue;
      for (int c = 0; c < q; ++c) {
        A[r * q + c] -= factor * A[col * q + c];
      }
      rhs[r] -= factor * rhs[col];
    }
  }
  return true;
}

// Weighted local linear regression of target on the numeric predictors that are
// observed for the receiver, fit on the selected donor neighborhood. Falls back
// to the caller's weighted-mean estimate (via ok = false) when the local design
// is degenerate (too few usable donors relative to predictors).
static bool missknn_local_regression(
    const RMatrix<double>& numeric_scaled,
    const RMatrix<double>& numeric_raw,
    const std::vector<int>& numeric_pos,
    const std::vector<int>& col_types,
    int target_col,
    int target_pos,
    int receiver_row,
    const std::vector<int>& donor_local_idx,
    const std::vector<int>& donor_rows,
    const std::vector<double>& donor_weights,
    double ridge,
    double& out_value,
    double& out_resid_sd) {
  const int p = static_cast<int>(col_types.size());
  std::vector<int> pred_pos;
  pred_pos.reserve(p);
  for (int col = 0; col < p; ++col) {
    if (col == target_col || col_types[col] != 1) {
      continue;
    }
    const int pos = numeric_pos[col] - 1;
    if (pos < 0) continue;
    const double rv = numeric_scaled(receiver_row - 1, pos);
    if (!R_IsNA(rv)) {
      pred_pos.push_back(pos);
    }
  }

  const int q = static_cast<int>(pred_pos.size()) + 1;
  const int n_don = static_cast<int>(donor_local_idx.size());
  if (n_don < q + 2) {
    return false;
  }

  std::vector<double> X(static_cast<size_t>(n_don) * q);
  std::vector<double> y(n_don);
  std::vector<double> w(n_don);
  int used = 0;
  for (int di = 0; di < n_don; ++di) {
    const int donor_row = donor_rows[donor_local_idx[di]] - 1;
    const double yv = numeric_raw(donor_row, target_pos);
    if (R_IsNA(yv)) continue;
    bool complete = true;
    std::vector<double> row(q);
    row[0] = 1.0;
    for (size_t j = 0; j < pred_pos.size(); ++j) {
      const double xv = numeric_scaled(donor_row, pred_pos[j]);
      if (R_IsNA(xv)) {
        complete = false;
        break;
      }
      row[j + 1] = xv;
    }
    if (!complete) continue;
    for (int c = 0; c < q; ++c) {
      X[static_cast<size_t>(used) * q + c] = row[c];
    }
    y[used] = yv;
    w[used] = donor_weights[di];
    ++used;
  }

  if (used < q + 2) {
    return false;
  }

  std::vector<double> XtWX(static_cast<size_t>(q) * q, 0.0);
  std::vector<double> XtWy(q, 0.0);
  for (int i = 0; i < used; ++i) {
    const double wi = w[i];
    for (int a = 0; a < q; ++a) {
      const double xa = X[static_cast<size_t>(i) * q + a];
      XtWy[a] += wi * xa * y[i];
      for (int b = 0; b < q; ++b) {
        XtWX[static_cast<size_t>(a) * q + b] += wi * xa * X[static_cast<size_t>(i) * q + b];
      }
    }
  }
  for (int a = 1; a < q; ++a) {
    XtWX[static_cast<size_t>(a) * q + a] += ridge;
  }

  std::vector<double> beta = XtWy;
  if (!solve_linear_system(XtWX, beta, q)) {
    return false;
  }

  double resid_ss = 0.0;
  double wsum = 0.0;
  for (int i = 0; i < used; ++i) {
    double fitted = 0.0;
    for (int c = 0; c < q; ++c) {
      fitted += beta[c] * X[static_cast<size_t>(i) * q + c];
    }
    const double e = y[i] - fitted;
    resid_ss += w[i] * e * e;
    wsum += w[i];
  }
  out_resid_sd = wsum > 0.0 ? std::sqrt(std::max(resid_ss / wsum, 0.0)) : 0.0;

  double pred = beta[0];
  for (size_t j = 0; j < pred_pos.size(); ++j) {
    pred += beta[j + 1] * numeric_scaled(receiver_row - 1, pred_pos[j]);
  }

  // A linear fit can extrapolate outside the observed donor range in a way a
  // tree-based model structurally cannot -- e.g. predicting a negative value
  // for a variable (like a follow-up duration) that is never negative in the
  // data. Clipping to the observed range of the donors actually used in this
  // fit removes that failure mode without changing the estimate whenever the
  // linear fit already falls inside the data's support, which is the common
  // case.
  double y_min = y[0], y_max = y[0];
  for (int i = 1; i < used; ++i) {
    if (y[i] < y_min) y_min = y[i];
    if (y[i] > y_max) y_max = y[i];
  }
  if (pred < y_min) pred = y_min;
  if (pred > y_max) pred = y_max;

  out_value = pred;
  return true;
}

// Imputes a single receiver's numeric value; shared by the serial (stochastic)
// loop and the parallel (deterministic) worker below.
static double missknn_impute_numeric_one(
    int receiver,
    const RMatrix<double>& numeric_scaled,
    const RMatrix<double>& numeric_raw,
    const RMatrix<int>& categorical_codes,
    const std::vector<int>& col_types,
    const std::vector<int>& numeric_pos,
    const std::vector<int>& categorical_pos,
    const std::vector<int>& donor_rows,
    int target_col,
    int target_pos,
    const std::vector<double>& weights,
    int k,
    const std::string& aggregation,
    const std::string& estimator,
    double epsilon,
    double ridge,
    int reg_neighbors,
    std::vector<double>& distances_buf,
    std::vector<int>& topk_buf,
    std::vector<int>& topk_reg_buf) {
  const int n_don = static_cast<int>(donor_rows.size());
  for (int j = 0; j < n_don; ++j) {
    distances_buf[j] = missknn_distance(receiver - 1, donor_rows[j] - 1, target_col - 1,
                                         numeric_scaled, categorical_codes, col_types,
                                         numeric_pos, categorical_pos, weights);
  }
  missknn_topk(distances_buf, donor_rows, k, topk_buf);

  if (topk_buf.empty()) {
    return missknn_numeric_fallback(numeric_raw, donor_rows, target_pos);
  }

  std::vector<double> sel_values;
  std::vector<double> sel_weights;
  sel_values.reserve(topk_buf.size());
  sel_weights.reserve(topk_buf.size());
  for (int idx : topk_buf) {
    const int donor = donor_rows[idx] - 1;
    const double d = distances_buf[idx];
    const double v = numeric_raw(donor, target_pos);
    if (R_IsNA(v)) continue;
    const double w = aggregation == "uniform" ? 1.0 : 1.0 / (d + epsilon);
    sel_values.push_back(v);
    sel_weights.push_back(w);
  }

  if (sel_values.empty()) {
    return missknn_numeric_fallback(numeric_raw, donor_rows, target_pos);
  }

  double mean_fallback_sum_w = 0.0, mean_fallback_sum_v = 0.0;
  for (size_t i = 0; i < sel_values.size(); ++i) {
    mean_fallback_sum_w += sel_weights[i];
    mean_fallback_sum_v += sel_weights[i] * sel_values[i];
  }
  const double mean_estimate = mean_fallback_sum_w > 0.0
    ? mean_fallback_sum_v / mean_fallback_sum_w
    : missknn_numeric_fallback(numeric_raw, donor_rows, target_pos);

  double point_estimate = mean_estimate;
  if (estimator == "regression") {
    missknn_topk(distances_buf, donor_rows, reg_neighbors, topk_reg_buf);
    // Uniform weights within the neighborhood: inverse-distance weights are
    // too peaked here (a near-duplicate donor dominates via leverage and the
    // fit can extrapolate wildly, unlike an averaged estimate). Locality is
    // already provided by restricting to the k-nearest neighborhood.
    std::vector<double> reg_weights(topk_reg_buf.size(), 1.0);
    double reg_value = 0.0;
    double reg_resid_sd = 0.0;
    const bool ok = missknn_local_regression(
        numeric_scaled, numeric_raw, numeric_pos, col_types,
        target_col - 1, target_pos, receiver, topk_reg_buf, donor_rows, reg_weights,
        ridge, reg_value, reg_resid_sd);
    if (ok && R_finite(reg_value)) {
      point_estimate = reg_value;
    }
  }

  return point_estimate;
}

// Parallel worker for the deterministic (non-stochastic) numeric receiver
// loop. Each receiver's distance search, top-k selection and prediction are
// fully independent of every other receiver, so this is embarrassingly
// parallel; the stochastic path (m > 1) is kept serial because it draws from
// R's RNG, which is not thread-safe.
struct ImputeNumericWorker : public Worker {
  const RMatrix<double> numeric_scaled;
  const RMatrix<double> numeric_raw;
  const RMatrix<int> categorical_codes;
  const std::vector<int>& col_types;
  const std::vector<int>& numeric_pos;
  const std::vector<int>& categorical_pos;
  const std::vector<int>& donor_rows;
  const std::vector<int>& receiver_rows;
  const int target_col;
  const int target_pos;
  const std::vector<double>& weights;
  const int k;
  const std::string& aggregation;
  const std::string& estimator;
  const double epsilon;
  const double ridge;
  const int reg_neighbors;
  RVector<double> out;

  ImputeNumericWorker(const NumericMatrix& numeric_scaled_,
                       const NumericMatrix& numeric_raw_,
                       const IntegerMatrix& categorical_codes_,
                       const std::vector<int>& col_types_,
                       const std::vector<int>& numeric_pos_,
                       const std::vector<int>& categorical_pos_,
                       const std::vector<int>& donor_rows_,
                       const std::vector<int>& receiver_rows_,
                       int target_col_,
                       int target_pos_,
                       const std::vector<double>& weights_,
                       int k_,
                       const std::string& aggregation_,
                       const std::string& estimator_,
                       double epsilon_,
                       double ridge_,
                       int reg_neighbors_,
                       NumericVector out_)
    : numeric_scaled(numeric_scaled_), numeric_raw(numeric_raw_), categorical_codes(categorical_codes_),
      col_types(col_types_), numeric_pos(numeric_pos_), categorical_pos(categorical_pos_),
      donor_rows(donor_rows_), receiver_rows(receiver_rows_), target_col(target_col_), target_pos(target_pos_),
      weights(weights_), k(k_), aggregation(aggregation_), estimator(estimator_),
      epsilon(epsilon_), ridge(ridge_), reg_neighbors(reg_neighbors_), out(out_) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> distances(donor_rows.size());
    std::vector<int> topk;
    std::vector<int> topk_reg;
    for (std::size_t r = begin; r < end; ++r) {
      out[r] = missknn_impute_numeric_one(
          receiver_rows[r], numeric_scaled, numeric_raw, categorical_codes,
          col_types, numeric_pos, categorical_pos, donor_rows, target_col, target_pos,
          weights, k, aggregation, estimator, epsilon, ridge, reg_neighbors,
          distances, topk, topk_reg);
    }
  }
};

// [[Rcpp::export]]
NumericVector cpp_missknn_impute_numeric_column(
    const NumericMatrix& numeric_scaled,
    const NumericMatrix& numeric_raw,
    const IntegerMatrix& categorical_codes,
    const IntegerVector& col_types,
    const IntegerVector& numeric_pos,
    const IntegerVector& categorical_pos,
    const IntegerVector& donor_rows,
    const IntegerVector& receiver_rows,
    const int target_col,
    const NumericVector& weights,
    const int k,
    const std::string& aggregation,
    const std::string& estimator,
    const bool stochastic,
    const double epsilon,
    const double ridge,
    const int reg_neighbors,
    const bool display_progress = false) {
  const int n_recv = receiver_rows.size();
  const int n_don = donor_rows.size();
  const int target_pos = numeric_pos[target_col - 1] - 1;
  NumericVector out(n_recv);

  if (n_don == 0) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }

  const std::vector<int> col_types_v(col_types.begin(), col_types.end());
  const std::vector<int> numeric_pos_v(numeric_pos.begin(), numeric_pos.end());
  const std::vector<int> categorical_pos_v(categorical_pos.begin(), categorical_pos.end());
  const std::vector<int> donor_rows_v(donor_rows.begin(), donor_rows.end());
  const std::vector<int> receiver_rows_v(receiver_rows.begin(), receiver_rows.end());
  const std::vector<double> weights_v(weights.begin(), weights.end());

  if (!stochastic) {
    ImputeNumericWorker worker(numeric_scaled, numeric_raw, categorical_codes,
                                col_types_v, numeric_pos_v, categorical_pos_v,
                                donor_rows_v, receiver_rows_v, target_col, target_pos,
                                weights_v, k, aggregation, estimator, epsilon, ridge,
                                reg_neighbors, out);
    missknn_run_with_optional_progress(worker, n_recv, display_progress);
    return out;
  }

  const RMatrix<double> ns(numeric_scaled);
  const RMatrix<double> nr(numeric_raw);
  const RMatrix<int> cc(categorical_codes);
  std::vector<double> distances(n_don);
  std::vector<int> topk;
  std::vector<int> topk_reg;

  Progress p_stoch(static_cast<unsigned long>(n_recv), display_progress);
  for (int r = 0; r < n_recv; ++r) {
    if (Progress::check_abort()) {
      Rcpp::stop("missknn: interrupted by user.");
    }
    p_stoch.increment();
    const int receiver = receiver_rows_v[r];
    for (int j = 0; j < n_don; ++j) {
      distances[j] = missknn_distance(receiver - 1, donor_rows_v[j] - 1, target_col - 1,
                                       ns, cc, col_types_v,
                                       numeric_pos_v, categorical_pos_v, weights_v);
    }
    missknn_topk(distances, donor_rows_v, k, topk);

    if (topk.empty()) {
      out[r] = missknn_numeric_fallback(nr, donor_rows_v, target_pos);
      continue;
    }

    std::vector<double> sel_values;
    std::vector<double> sel_weights;
    sel_values.reserve(topk.size());
    sel_weights.reserve(topk.size());
    for (int idx : topk) {
      const int donor = donor_rows_v[idx] - 1;
      const double d = distances[idx];
      const double v = nr(donor, target_pos);
      if (R_IsNA(v)) continue;
      const double w = aggregation == "uniform" ? 1.0 : 1.0 / (d + epsilon);
      sel_values.push_back(v);
      sel_weights.push_back(w);
    }

    if (sel_values.empty()) {
      out[r] = missknn_numeric_fallback(nr, donor_rows_v, target_pos);
      continue;
    }

    double mean_fallback_sum_w = 0.0, mean_fallback_sum_v = 0.0;
    for (size_t i = 0; i < sel_values.size(); ++i) {
      mean_fallback_sum_w += sel_weights[i];
      mean_fallback_sum_v += sel_weights[i] * sel_values[i];
    }
    const double mean_estimate = mean_fallback_sum_w > 0.0
      ? mean_fallback_sum_v / mean_fallback_sum_w
      : missknn_numeric_fallback(nr, donor_rows_v, target_pos);

    double point_estimate = mean_estimate;
    double resid_sd = 0.0;
    if (estimator == "regression") {
      missknn_topk(distances, donor_rows_v, reg_neighbors, topk_reg);
      std::vector<double> reg_weights(topk_reg.size(), 1.0);
      double reg_value = 0.0;
      double reg_resid_sd = 0.0;
      const bool ok = missknn_local_regression(
          ns, nr, numeric_pos_v, col_types_v,
          target_col - 1, target_pos, receiver, topk_reg, donor_rows_v, reg_weights,
          ridge, reg_value, reg_resid_sd);
      if (ok && R_finite(reg_value)) {
        point_estimate = reg_value;
        resid_sd = reg_resid_sd;
      }
    }

    if ((estimator == "regression") && resid_sd > 0.0) {
      out[r] = point_estimate + R::rnorm(0.0, resid_sd);
    } else {
      double sum_w = 0.0;
      for (double w : sel_weights) sum_w += w;
      if (sum_w <= 0.0) {
        const int idx = static_cast<int>(std::floor(R::runif(0.0, 1.0) * sel_values.size()));
        out[r] = sel_values[std::min(idx, static_cast<int>(sel_values.size()) - 1)];
      } else {
        const double u = R::runif(0.0, sum_w);
        double acc = 0.0;
        double picked = sel_values.back();
        for (size_t i = 0; i < sel_values.size(); ++i) {
          acc += sel_weights[i];
          if (u <= acc) { picked = sel_values[i]; break; }
        }
        out[r] = picked;
      }
    }
  }

  return out;
}

// Imputes a single receiver's categorical code (majority-vote branch only;
// the stochastic sampling branch stays in the serial caller below since it
// draws from R's RNG).
static int missknn_impute_categorical_one(
    int receiver,
    const RMatrix<double>& numeric_scaled,
    const RMatrix<int>& categorical_codes,
    const std::vector<int>& col_types,
    const std::vector<int>& numeric_pos,
    const std::vector<int>& categorical_pos,
    const std::vector<int>& donor_rows,
    int target_col,
    int target_pos,
    const std::vector<double>& weights,
    int k,
    const std::string& aggregation,
    double epsilon,
    std::vector<double>& distances_buf,
    std::vector<int>& topk_buf) {
  const int n_don = static_cast<int>(donor_rows.size());
  for (int j = 0; j < n_don; ++j) {
    distances_buf[j] = missknn_distance(receiver - 1, donor_rows[j] - 1, target_col - 1,
                                         numeric_scaled, categorical_codes, col_types,
                                         numeric_pos, categorical_pos, weights);
  }
  missknn_topk(distances_buf, donor_rows, k, topk_buf);

  if (topk_buf.empty()) {
    return missknn_categorical_fallback(categorical_codes, donor_rows, target_pos);
  }

  std::vector<int> sel_codes;
  std::vector<double> sel_weights;
  sel_codes.reserve(topk_buf.size());
  sel_weights.reserve(topk_buf.size());
  for (int idx : topk_buf) {
    const int donor = donor_rows[idx] - 1;
    const double d = distances_buf[idx];
    const int code = categorical_codes(donor, target_pos);
    if (code == NA_INTEGER) continue;
    const double w = aggregation == "uniform" ? 1.0 : 1.0 / (d + epsilon);
    sel_codes.push_back(code);
    sel_weights.push_back(w);
  }

  if (sel_codes.empty()) {
    return missknn_categorical_fallback(categorical_codes, donor_rows, target_pos);
  }

  std::map<int, double> scores;
  for (size_t i = 0; i < sel_codes.size(); ++i) {
    scores[sel_codes[i]] += sel_weights[i];
  }
  int best_code = NA_INTEGER;
  double best_score = R_NegInf;
  for (const auto& kv : scores) {
    if (kv.second > best_score || (kv.second == best_score && kv.first < best_code)) {
      best_code = kv.first;
      best_score = kv.second;
    }
  }
  return best_code;
}

struct ImputeCategoricalWorker : public Worker {
  const RMatrix<double> numeric_scaled;
  const RMatrix<int> categorical_codes;
  const std::vector<int>& col_types;
  const std::vector<int>& numeric_pos;
  const std::vector<int>& categorical_pos;
  const std::vector<int>& donor_rows;
  const std::vector<int>& receiver_rows;
  const int target_col;
  const int target_pos;
  const std::vector<double>& weights;
  const int k;
  const std::string& aggregation;
  const double epsilon;
  RVector<int> out;

  ImputeCategoricalWorker(const NumericMatrix& numeric_scaled_,
                           const IntegerMatrix& categorical_codes_,
                           const std::vector<int>& col_types_,
                           const std::vector<int>& numeric_pos_,
                           const std::vector<int>& categorical_pos_,
                           const std::vector<int>& donor_rows_,
                           const std::vector<int>& receiver_rows_,
                           int target_col_,
                           int target_pos_,
                           const std::vector<double>& weights_,
                           int k_,
                           const std::string& aggregation_,
                           double epsilon_,
                           IntegerVector out_)
    : numeric_scaled(numeric_scaled_), categorical_codes(categorical_codes_),
      col_types(col_types_), numeric_pos(numeric_pos_), categorical_pos(categorical_pos_),
      donor_rows(donor_rows_), receiver_rows(receiver_rows_), target_col(target_col_), target_pos(target_pos_),
      weights(weights_), k(k_), aggregation(aggregation_), epsilon(epsilon_), out(out_) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> distances(donor_rows.size());
    std::vector<int> topk;
    for (std::size_t r = begin; r < end; ++r) {
      out[r] = missknn_impute_categorical_one(
          receiver_rows[r], numeric_scaled, categorical_codes, col_types,
          numeric_pos, categorical_pos, donor_rows, target_col, target_pos,
          weights, k, aggregation, epsilon, distances, topk);
    }
  }
};

// [[Rcpp::export]]
IntegerVector cpp_missknn_impute_categorical_column(
    const NumericMatrix& numeric_scaled,
    const IntegerMatrix& categorical_codes,
    const IntegerVector& col_types,
    const IntegerVector& numeric_pos,
    const IntegerVector& categorical_pos,
    const IntegerVector& donor_rows,
    const IntegerVector& receiver_rows,
    const int target_col,
    const NumericVector& weights,
    const int k,
    const std::string& aggregation,
    const bool stochastic,
    const double epsilon,
    const bool display_progress = false) {
  const int n_recv = receiver_rows.size();
  const int n_don = donor_rows.size();
  const int target_pos = categorical_pos[target_col - 1] - 1;
  IntegerVector out(n_recv);

  if (n_don == 0) {
    std::fill(out.begin(), out.end(), NA_INTEGER);
    return out;
  }

  const std::vector<int> col_types_v(col_types.begin(), col_types.end());
  const std::vector<int> numeric_pos_v(numeric_pos.begin(), numeric_pos.end());
  const std::vector<int> categorical_pos_v(categorical_pos.begin(), categorical_pos.end());
  const std::vector<int> donor_rows_v(donor_rows.begin(), donor_rows.end());
  const std::vector<int> receiver_rows_v(receiver_rows.begin(), receiver_rows.end());
  const std::vector<double> weights_v(weights.begin(), weights.end());

  if (!stochastic) {
    ImputeCategoricalWorker worker(numeric_scaled, categorical_codes, col_types_v,
                                    numeric_pos_v, categorical_pos_v, donor_rows_v,
                                    receiver_rows_v, target_col, target_pos, weights_v,
                                    k, aggregation, epsilon, out);
    missknn_run_with_optional_progress(worker, n_recv, display_progress);
    return out;
  }

  const RMatrix<double> ns(numeric_scaled);
  const RMatrix<int> cc(categorical_codes);
  std::vector<double> distances(n_don);
  std::vector<int> topk;

  Progress p_stoch(static_cast<unsigned long>(n_recv), display_progress);
  for (int r = 0; r < n_recv; ++r) {
    if (Progress::check_abort()) {
      Rcpp::stop("missknn: interrupted by user.");
    }
    p_stoch.increment();
    const int receiver = receiver_rows_v[r];
    for (int j = 0; j < n_don; ++j) {
      distances[j] = missknn_distance(receiver - 1, donor_rows_v[j] - 1, target_col - 1,
                                       ns, cc, col_types_v,
                                       numeric_pos_v, categorical_pos_v, weights_v);
    }
    missknn_topk(distances, donor_rows_v, k, topk);

    if (topk.empty()) {
      out[r] = missknn_categorical_fallback(cc, donor_rows_v, target_pos);
      continue;
    }

    std::vector<int> sel_codes;
    std::vector<double> sel_weights;
    sel_codes.reserve(topk.size());
    sel_weights.reserve(topk.size());
    for (int idx : topk) {
      const int donor = donor_rows_v[idx] - 1;
      const double d = distances[idx];
      const int code = cc(donor, target_pos);
      if (code == NA_INTEGER) continue;
      const double w = aggregation == "uniform" ? 1.0 : 1.0 / (d + epsilon);
      sel_codes.push_back(code);
      sel_weights.push_back(w);
    }

    if (sel_codes.empty()) {
      out[r] = missknn_categorical_fallback(cc, donor_rows_v, target_pos);
      continue;
    }

    double sum_w = 0.0;
    for (double w : sel_weights) sum_w += w;
    if (sum_w <= 0.0) {
      const int idx = static_cast<int>(std::floor(R::runif(0.0, 1.0) * sel_codes.size()));
      out[r] = sel_codes[std::min(idx, static_cast<int>(sel_codes.size()) - 1)];
    } else {
      const double u = R::runif(0.0, sum_w);
      double acc = 0.0;
      int picked = sel_codes.back();
      for (size_t i = 0; i < sel_codes.size(); ++i) {
        acc += sel_weights[i];
        if (u <= acc) { picked = sel_codes[i]; break; }
      }
      out[r] = picked;
    }
  }

  return out;
}

// Per-column result of the numeric tuning search (see cpp_missknn_tune_numeric).
struct NumericTuneResult {
  int k = 1;
  int estimator_code = 0; // 0 = "mean", 1 = "regression"
  bool is_global = true;
};

// Tunes k (and mean vs. local-regression estimator) for one numeric target
// column: caches the hold x train masked-distance matrix and reuses it across
// all candidate k values and the estimator check, avoiding both repeated
// distance computation and per-configuration call overhead. Columns are
// fully independent of one another, which is what lets the caller run this
// once per column across threads.
static NumericTuneResult missknn_tune_numeric_one_column(
    const RMatrix<double>& numeric_scaled,
    const RMatrix<double>& numeric_raw,
    const RMatrix<int>& categorical_codes,
    const std::vector<int>& col_types,
    const std::vector<int>& numeric_pos,
    const std::vector<int>& categorical_pos,
    const std::vector<double>& weights,
    double epsilon,
    double ridge,
    int target_col,
    const std::vector<int>& train_idx,
    const std::vector<int>& hold_idx,
    const std::vector<int>& k_grid) {
  const int target_pos = numeric_pos[target_col - 1] - 1;
  const int n_train = static_cast<int>(train_idx.size());
  const int n_hold = static_cast<int>(hold_idx.size());

  std::vector<std::vector<double>> dist_mat(n_hold, std::vector<double>(n_train));
  for (int h = 0; h < n_hold; ++h) {
    for (int j = 0; j < n_train; ++j) {
      dist_mat[h][j] = missknn_distance(hold_idx[h] - 1, train_idx[j] - 1, target_col - 1,
                                         numeric_scaled, categorical_codes, col_types,
                                         numeric_pos, categorical_pos, weights);
    }
  }

  std::vector<double> truth(n_hold);
  for (int h = 0; h < n_hold; ++h) {
    truth[h] = numeric_raw(hold_idx[h] - 1, target_pos);
  }

  // Rank donors once per holdout point (full order); every candidate k and
  // the regression neighborhood just take a prefix of this ranking, so the
  // O(n_train) selection cost is paid once instead of once per grid point.
  std::vector<std::vector<int>> ranked(n_hold);
  for (int h = 0; h < n_hold; ++h) {
    missknn_topk(dist_mat[h], train_idx, n_train, ranked[h]);
  }

  // Baseline candidate: the plain (unweighted) global mean over all donors.
  // When a column carries no locally exploitable signal (e.g. independent
  // noise), this beats any distance-weighted neighborhood average, since
  // inverse-distance weights add avoidable variance without adding bias
  // correction. Columns that keep this baseline skip the O(receivers x
  // donors) KNN machinery entirely at imputation time.
  double global_sum = 0.0;
  int global_n = 0;
  for (int j = 0; j < n_train; ++j) {
    const double v = numeric_raw(train_idx[j] - 1, target_pos);
    if (!R_IsNA(v)) {
      global_sum += v;
      ++global_n;
    }
  }
  const double global_mean = global_n > 0 ? global_sum / global_n : NA_REAL;
  double best_err = R_PosInf;
  if (global_n > 0) {
    double sse = 0.0;
    int cnt = 0;
    for (int h = 0; h < n_hold; ++h) {
      if (R_IsNA(truth[h])) continue;
      const double e = global_mean - truth[h];
      sse += e * e;
      ++cnt;
    }
    if (cnt > 0) {
      best_err = sse / cnt;
    }
  }
  int chosen_k = k_grid.size() ? k_grid[0] : 1;
  int chosen_est_code = 0;
  bool chosen_global = true;

  for (size_t ki = 0; ki < k_grid.size(); ++ki) {
    const int kk = k_grid[ki];
    double sse = 0.0;
    int cnt = 0;
    for (int h = 0; h < n_hold; ++h) {
      if (R_IsNA(truth[h])) continue;
      const int use_k = std::min(kk, static_cast<int>(ranked[h].size()));
      if (use_k <= 0) continue;
      double sw = 0.0, sv = 0.0;
      for (int ii = 0; ii < use_k; ++ii) {
        const int idx = ranked[h][ii];
        const double v = numeric_raw(train_idx[idx] - 1, target_pos);
        if (R_IsNA(v)) continue;
        const double w = 1.0 / (dist_mat[h][idx] + epsilon);
        sw += w;
        sv += w * v;
      }
      if (sw <= 0.0) continue;
      const double e = (sv / sw) - truth[h];
      sse += e * e;
      ++cnt;
    }
    if (cnt > 0) {
      const double mse = sse / cnt;
      // Require a real relative improvement over the current best (initially
      // the global baseline) before switching, not just a nominally lower
      // holdout error: with a finite holdout sample, a local k can look
      // marginally better than the global mean purely by chance even when a
      // column carries no exploitable local signal, which would otherwise
      // send it down the slower, noisier KNN path for no real accuracy gain.
      if (mse < best_err * (1.0 - kTuneMargin)) {
        best_err = mse;
        chosen_k = kk;
        chosen_global = false;
      }
    }
  }

  // Exactly two regression bandwidth candidates -- narrow (a small multiple
  // of the winning mean-estimator k) and wide (capped well below the full
  // donor pool, since the variance-reduction benefit of more donors
  // saturates long before that and an unbounded wide fit has no O(1) fast
  // path if a false positive sends every receiver through it). Kept to two
  // candidates deliberately: picking the minimum error over several noisy
  // candidates is itself biased low even under pure noise (the same effect
  // that makes "best of several dice rolls" beat a fair average), so more
  // candidates raises the false-positive rate before the margin check even
  // runs; kRegressionMargin compensates for the bias these two still carry.
  const double pre_reg_best_err = best_err;
  int best_reg_k = -1;
  double best_reg_err = R_PosInf;

  const int reg_base_k = chosen_global ? (k_grid.size() ? k_grid[k_grid.size() - 1] : 30) : chosen_k;
  const int reg_k_candidates[2] = {
    std::max(4 * reg_base_k, 30),
    std::min(n_train, 300)
  };

  for (size_t ci = 0; ci < 2; ++ci) {
    const int reg_k = reg_k_candidates[ci];
    double sse = 0.0;
    int cnt = 0;
    for (int h = 0; h < n_hold; ++h) {
      if (R_IsNA(truth[h])) continue;
      const int use_k = std::min(reg_k, static_cast<int>(ranked[h].size()));
      if (use_k <= 0) continue;
      std::vector<int> topk_reg(ranked[h].begin(), ranked[h].begin() + use_k);
      std::vector<double> rw(topk_reg.size(), 1.0);
      double val = 0.0, resid_sd = 0.0;
      const bool ok = missknn_local_regression(
          numeric_scaled, numeric_raw, numeric_pos, col_types,
          target_col - 1, target_pos, hold_idx[h], topk_reg, train_idx, rw,
          ridge, val, resid_sd);
      if (!ok || !R_finite(val)) continue;
      const double e = val - truth[h];
      sse += e * e;
      ++cnt;
    }
    if (cnt == 0) continue;
    const double mse = sse / cnt;
    if (mse < best_reg_err) {
      best_reg_err = mse;
      best_reg_k = reg_k;
    }
  }

  if (best_reg_k > 0 && best_reg_err < pre_reg_best_err * (1.0 - kRegressionMargin)) {
    chosen_est_code = 1;
    chosen_k = best_reg_k;
    chosen_global = false;
  }

  NumericTuneResult result;
  result.k = chosen_k;
  result.estimator_code = chosen_est_code;
  result.is_global = chosen_global;
  return result;
}

// Parallel worker over target columns: each column's tuning search reads
// the shared matrices but writes only to its own slot t of the output
// arrays, so columns can run concurrently with no shared mutable state.
struct TuneNumericWorker : public Worker {
  const RMatrix<double> numeric_scaled;
  const RMatrix<double> numeric_raw;
  const RMatrix<int> categorical_codes;
  const std::vector<int>& col_types;
  const std::vector<int>& numeric_pos;
  const std::vector<int>& categorical_pos;
  const std::vector<double>& weights;
  const double epsilon;
  const double ridge;
  const std::vector<int>& target_cols;
  const std::vector<std::vector<int>>& train_list;
  const std::vector<std::vector<int>>& hold_list;
  const std::vector<std::vector<int>>& k_grid_list;
  std::vector<int>& out_k;
  std::vector<int>& out_est;
  std::vector<int>& out_global;

  TuneNumericWorker(const NumericMatrix& numeric_scaled_,
                     const NumericMatrix& numeric_raw_,
                     const IntegerMatrix& categorical_codes_,
                     const std::vector<int>& col_types_,
                     const std::vector<int>& numeric_pos_,
                     const std::vector<int>& categorical_pos_,
                     const std::vector<double>& weights_,
                     double epsilon_,
                     double ridge_,
                     const std::vector<int>& target_cols_,
                     const std::vector<std::vector<int>>& train_list_,
                     const std::vector<std::vector<int>>& hold_list_,
                     const std::vector<std::vector<int>>& k_grid_list_,
                     std::vector<int>& out_k_,
                     std::vector<int>& out_est_,
                     std::vector<int>& out_global_)
    : numeric_scaled(numeric_scaled_), numeric_raw(numeric_raw_), categorical_codes(categorical_codes_),
      col_types(col_types_), numeric_pos(numeric_pos_), categorical_pos(categorical_pos_), weights(weights_),
      epsilon(epsilon_), ridge(ridge_), target_cols(target_cols_), train_list(train_list_),
      hold_list(hold_list_), k_grid_list(k_grid_list_), out_k(out_k_), out_est(out_est_), out_global(out_global_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t t = begin; t < end; ++t) {
      const NumericTuneResult res = missknn_tune_numeric_one_column(
          numeric_scaled, numeric_raw, categorical_codes, col_types, numeric_pos, categorical_pos,
          weights, epsilon, ridge, target_cols[t], train_list[t], hold_list[t], k_grid_list[t]);
      out_k[t] = res.k;
      out_est[t] = res.estimator_code;
      out_global[t] = res.is_global ? 1 : 0;
    }
  }
};

// [[Rcpp::export]]
List cpp_missknn_tune_numeric(
    const NumericMatrix& numeric_scaled,
    const NumericMatrix& numeric_raw,
    const IntegerMatrix& categorical_codes,
    const IntegerVector& col_types,
    const IntegerVector& numeric_pos,
    const IntegerVector& categorical_pos,
    const NumericVector& weights,
    const double epsilon,
    const double ridge,
    const IntegerVector& target_cols,
    const List& train_list,
    const List& hold_list,
    const List& k_grid_list,
    const bool display_progress = false) {
  const int ntar = target_cols.size();

  const std::vector<int> col_types_v(col_types.begin(), col_types.end());
  const std::vector<int> numeric_pos_v(numeric_pos.begin(), numeric_pos.end());
  const std::vector<int> categorical_pos_v(categorical_pos.begin(), categorical_pos.end());
  const std::vector<double> weights_v(weights.begin(), weights.end());
  const std::vector<int> target_cols_v(target_cols.begin(), target_cols.end());

  std::vector<std::vector<int>> train_v(ntar), hold_v(ntar), k_grid_v(ntar);
  for (int t = 0; t < ntar; ++t) {
    const IntegerVector tr = train_list[t];
    const IntegerVector ho = hold_list[t];
    const IntegerVector kg = k_grid_list[t];
    train_v[t].assign(tr.begin(), tr.end());
    hold_v[t].assign(ho.begin(), ho.end());
    k_grid_v[t].assign(kg.begin(), kg.end());
  }

  std::vector<int> out_k(ntar), out_est(ntar), out_global(ntar);

  TuneNumericWorker worker(numeric_scaled, numeric_raw, categorical_codes, col_types_v,
                            numeric_pos_v, categorical_pos_v, weights_v, epsilon, ridge,
                            target_cols_v, train_v, hold_v, k_grid_v, out_k, out_est, out_global);
  missknn_run_with_optional_progress(worker, ntar, display_progress);

  IntegerVector best_k(ntar);
  CharacterVector best_est(ntar);
  LogicalVector is_global(ntar);
  for (int t = 0; t < ntar; ++t) {
    best_k[t] = out_k[t];
    best_est[t] = out_est[t] == 1 ? "regression" : "mean";
    is_global[t] = out_global[t] == 1;
  }

  return List::create(Named("k") = best_k, Named("estimator") = best_est, Named("is_global") = is_global);
}

struct CategoricalTuneResult {
  int k = 1;
  bool is_global = true;
};

// Analogous per-column tuning for a categorical target column, selecting k
// that minimizes holdout misclassification rate.
static CategoricalTuneResult missknn_tune_categorical_one_column(
    const RMatrix<double>& numeric_scaled,
    const RMatrix<int>& categorical_codes,
    const std::vector<int>& col_types,
    const std::vector<int>& numeric_pos,
    const std::vector<int>& categorical_pos,
    const std::vector<double>& weights,
    double epsilon,
    int target_col,
    const std::vector<int>& train_idx,
    const std::vector<int>& hold_idx,
    const std::vector<int>& k_grid) {
  const int target_pos = categorical_pos[target_col - 1] - 1;
  const int n_train = static_cast<int>(train_idx.size());
  const int n_hold = static_cast<int>(hold_idx.size());

  std::vector<std::vector<double>> dist_mat(n_hold, std::vector<double>(n_train));
  for (int h = 0; h < n_hold; ++h) {
    for (int j = 0; j < n_train; ++j) {
      dist_mat[h][j] = missknn_distance(hold_idx[h] - 1, train_idx[j] - 1, target_col - 1,
                                         numeric_scaled, categorical_codes, col_types,
                                         numeric_pos, categorical_pos, weights);
    }
  }

  std::vector<int> truth(n_hold);
  for (int h = 0; h < n_hold; ++h) {
    truth[h] = categorical_codes(hold_idx[h] - 1, target_pos);
  }

  std::vector<std::vector<int>> ranked(n_hold);
  for (int h = 0; h < n_hold; ++h) {
    missknn_topk(dist_mat[h], train_idx, n_train, ranked[h]);
  }

  // Baseline candidate: the global mode over all donors.
  std::map<int, int> global_counts;
  for (int j = 0; j < n_train; ++j) {
    const int code = categorical_codes(train_idx[j] - 1, target_pos);
    if (code != NA_INTEGER) global_counts[code] += 1;
  }
  int global_mode = NA_INTEGER;
  int global_best_count = -1;
  for (const auto& kv : global_counts) {
    if (kv.second > global_best_count) {
      global_best_count = kv.second;
      global_mode = kv.first;
    }
  }
  double best_err = R_PosInf;
  if (global_mode != NA_INTEGER) {
    int wrong = 0, cnt = 0;
    for (int h = 0; h < n_hold; ++h) {
      if (truth[h] == NA_INTEGER) continue;
      ++cnt;
      if (global_mode != truth[h]) ++wrong;
    }
    if (cnt > 0) best_err = static_cast<double>(wrong) / cnt;
  }
  int chosen_k = k_grid.size() ? k_grid[0] : 1;
  bool chosen_global = true;

  auto eval_vote = [&](int kk) -> double {
    int wrong = 0;
    int cnt = 0;
    for (int h = 0; h < n_hold; ++h) {
      if (truth[h] == NA_INTEGER) continue;
      const int use_k = std::min(kk, static_cast<int>(ranked[h].size()));
      if (use_k <= 0) continue;
      std::map<int, double> scores;
      for (int ii = 0; ii < use_k; ++ii) {
        const int idx = ranked[h][ii];
        const int code = categorical_codes(train_idx[idx] - 1, target_pos);
        if (code == NA_INTEGER) continue;
        scores[code] += 1.0 / (dist_mat[h][idx] + epsilon);
      }
      if (scores.empty()) continue;
      int best_code = NA_INTEGER;
      double best_score = R_NegInf;
      for (const auto& kv : scores) {
        if (kv.second > best_score) {
          best_score = kv.second;
          best_code = kv.first;
        }
      }
      ++cnt;
      if (best_code != truth[h]) ++wrong;
    }
    return cnt > 0 ? static_cast<double>(wrong) / cnt : R_PosInf;
  };

  for (size_t ki = 0; ki < k_grid.size(); ++ki) {
    const int kk = k_grid[ki];
    const double err = eval_vote(kk);
    if (err < best_err * (1.0 - kTuneMargin)) {
      best_err = err;
      chosen_k = kk;
      chosen_global = false;
    }
  }

  // As with the numeric wide-regression candidate: a small local vote can
  // have needlessly high variance when the category is actually predicted
  // by a stable pattern across the whole donor pool rather than a locally
  // varying one. Evaluated against the pre-vote baseline directly (not
  // gated behind the narrow vote's own margin) for the same reason. The
  // chosen k is set far beyond any realistic donor pool (rather than to
  // this tuning split's own n_train) so it means "use every available
  // donor" at imputation time too, when the real donor pool -- capped at
  // donor_cap -- is typically larger than the tuning split.
  const int kWideSentinel = 1000000;
  const double err_wide = eval_vote(n_train);
  if (err_wide < best_err * (1.0 - kTuneMargin)) {
    best_err = err_wide;
    chosen_k = kWideSentinel;
    chosen_global = false;
  }

  CategoricalTuneResult result;
  result.k = chosen_k;
  result.is_global = chosen_global;
  return result;
}

struct TuneCategoricalWorker : public Worker {
  const RMatrix<double> numeric_scaled;
  const RMatrix<int> categorical_codes;
  const std::vector<int>& col_types;
  const std::vector<int>& numeric_pos;
  const std::vector<int>& categorical_pos;
  const std::vector<double>& weights;
  const double epsilon;
  const std::vector<int>& target_cols;
  const std::vector<std::vector<int>>& train_list;
  const std::vector<std::vector<int>>& hold_list;
  const std::vector<std::vector<int>>& k_grid_list;
  std::vector<int>& out_k;
  std::vector<int>& out_global;

  TuneCategoricalWorker(const NumericMatrix& numeric_scaled_,
                         const IntegerMatrix& categorical_codes_,
                         const std::vector<int>& col_types_,
                         const std::vector<int>& numeric_pos_,
                         const std::vector<int>& categorical_pos_,
                         const std::vector<double>& weights_,
                         double epsilon_,
                         const std::vector<int>& target_cols_,
                         const std::vector<std::vector<int>>& train_list_,
                         const std::vector<std::vector<int>>& hold_list_,
                         const std::vector<std::vector<int>>& k_grid_list_,
                         std::vector<int>& out_k_,
                         std::vector<int>& out_global_)
    : numeric_scaled(numeric_scaled_), categorical_codes(categorical_codes_), col_types(col_types_),
      numeric_pos(numeric_pos_), categorical_pos(categorical_pos_), weights(weights_), epsilon(epsilon_),
      target_cols(target_cols_), train_list(train_list_), hold_list(hold_list_), k_grid_list(k_grid_list_),
      out_k(out_k_), out_global(out_global_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t t = begin; t < end; ++t) {
      const CategoricalTuneResult res = missknn_tune_categorical_one_column(
          numeric_scaled, categorical_codes, col_types, numeric_pos, categorical_pos,
          weights, epsilon, target_cols[t], train_list[t], hold_list[t], k_grid_list[t]);
      out_k[t] = res.k;
      out_global[t] = res.is_global ? 1 : 0;
    }
  }
};

// [[Rcpp::export]]
List cpp_missknn_tune_categorical(
    const NumericMatrix& numeric_scaled,
    const IntegerMatrix& categorical_codes,
    const IntegerVector& col_types,
    const IntegerVector& numeric_pos,
    const IntegerVector& categorical_pos,
    const NumericVector& weights,
    const double epsilon,
    const IntegerVector& target_cols,
    const List& train_list,
    const List& hold_list,
    const List& k_grid_list,
    const bool display_progress = false) {
  const int ntar = target_cols.size();

  const std::vector<int> col_types_v(col_types.begin(), col_types.end());
  const std::vector<int> numeric_pos_v(numeric_pos.begin(), numeric_pos.end());
  const std::vector<int> categorical_pos_v(categorical_pos.begin(), categorical_pos.end());
  const std::vector<double> weights_v(weights.begin(), weights.end());
  const std::vector<int> target_cols_v(target_cols.begin(), target_cols.end());

  std::vector<std::vector<int>> train_v(ntar), hold_v(ntar), k_grid_v(ntar);
  for (int t = 0; t < ntar; ++t) {
    const IntegerVector tr = train_list[t];
    const IntegerVector ho = hold_list[t];
    const IntegerVector kg = k_grid_list[t];
    train_v[t].assign(tr.begin(), tr.end());
    hold_v[t].assign(ho.begin(), ho.end());
    k_grid_v[t].assign(kg.begin(), kg.end());
  }

  std::vector<int> out_k(ntar), out_global(ntar);

  TuneCategoricalWorker worker(numeric_scaled, categorical_codes, col_types_v, numeric_pos_v,
                                categorical_pos_v, weights_v, epsilon, target_cols_v,
                                train_v, hold_v, k_grid_v, out_k, out_global);
  missknn_run_with_optional_progress(worker, ntar, display_progress);

  IntegerVector best_k(ntar);
  LogicalVector is_global(ntar);
  for (int t = 0; t < ntar; ++t) {
    best_k[t] = out_k[t];
    is_global[t] = out_global[t] == 1;
  }

  return List::create(Named("k") = best_k, Named("is_global") = is_global);
}
