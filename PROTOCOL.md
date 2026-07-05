# missknn Protocol

## Purpose

`missknn` is a lightweight R package for KNN-based imputation of missing values in tabular data.

The design goal is not to replicate a full-featured imputation framework. The package should instead provide a fast, simple, and transparent engine for single imputation and multiple imputation with a minimal API.

## Positioning

`VIM::kNN()` is flexible and mature. `missknn` should be faster, smaller, more modern, and easier to use.

The core public claim for v1 is:

> `missknn` v1 will be a C++-backed, fast masked KNN imputation engine for R.

The core value proposition is:

- compact interface
- masked distance computation
- efficient donor selection
- reproducible single imputation
- multiple imputation support
- clean metadata around the completed dataset

## Non-Goals

`missknn` is not intended to be:

- a general-purpose hotdeck framework
- a random-forest-driven imputation suite
- a feature-heavy wrapper around many imputation methods
- a full MICE-style chained-equation framework

## Algorithm Name

Internal algorithm name:

**FM-KNN: Fast Masked k-Nearest Neighbor Imputation**

## Core Idea

Distances are computed only on jointly observed variables.

For a receiver row `i` and donor row `j`:

- only variables observed in both rows contribute to the distance
- missing values do not automatically penalize the pair
- the distance is normalized by the amount of shared information

This keeps comparisons fair and avoids over-penalizing rows with sparse but usable overlap.

## Method

This section defines the estimator in a form suitable for implementation and for a methods manuscript. The method has two modes:

- single imputation, which returns one completed dataset
- multiple imputation, which returns a collection of completed datasets generated from the same masked KNN mechanism

### Notation

Let the data be a tabular matrix or data frame with `n` rows and `p` variables:

```text
X = (x_{ij}) ∈ R^{n×p}
```

For each variable `j`, define the missingness indicator:

```text
R_{ij} = 1  if x_{ij} is observed
R_{ij} = 0  if x_{ij} is missing
```

Let `O_i = { j : R_{ij} = 1 }` be the set of observed variables for row `i`.

For a target variable `t` and a receiver row `i` with `R_{it} = 0`, define the donor set:

```text
D_t = { j : R_{jt} = 1 }
```

For any receiver row `i` and donor row `ℓ ∈ D_t`, define the common observed set excluding the target variable:

```text
S_{iℓt} = O_i ∩ O_ℓ ∩ {1, …, p} \ {t}
```

If `S_{iℓt} = ∅`, the pair is not comparable and is discarded.

### Variable Scaling

For numeric variables, let the scaled value be:

```text
z_{ij} = (x_{ij} - μ_j) / σ_j
```

with `μ_j` and `σ_j` estimated from observed values of variable `j`.

For categorical variables, use a finite encoding `c(x_{ij})` only for distance computation; the original labels are restored after imputation.

### Masked Distance Estimator

For a numeric-only path, define the masked squared distance between receiver `i` and donor `ℓ` for target `t` as:

```text
d^2_t(i, ℓ) =
  [ Σ_{j ∈ S_{iℓt}} w_j (z_{ij} - z_{ℓj})^2 ] /
  [ Σ_{j ∈ S_{iℓt}} w_j ]
```

and

```text
d_t(i, ℓ) = sqrt(d^2_t(i, ℓ)).
```

Here `w_j ≥ 0` is the variable weight and the denominator normalizes by the shared observed mass.

For mixed data, define a per-variable contribution `δ_j(i, ℓ)` by type:

```text
δ_j(i, ℓ) =
  (z_{ij} - z_{ℓj})^2                     if variable j is numeric
  1[x_{ij} ≠ x_{ℓj}]                       if variable j is categorical
```

Then the general masked distance is:

```text
d^2_t(i, ℓ) =
  [ Σ_{j ∈ S_{iℓt}} w_j δ_j(i, ℓ) ] /
  [ Σ_{j ∈ S_{iℓt}} w_j ]
```

### Neighbor Selection

For each receiver row `i` and target variable `t`, compute distances to all donors `ℓ ∈ D_t`.

Let `N_k(i, t)` be the index set of the `k` donors with smallest finite distance:

```text
N_k(i, t) = argmin_k { d_t(i, ℓ) : ℓ ∈ D_t, d_t(i, ℓ) < ∞ }
```

Ties are broken deterministically.

### Numeric Estimator

For a numeric target variable `t`, define donor weights:

```text
α_{iℓt} = 1 / (d_t(i, ℓ) + ε)
```

for `ℓ ∈ N_k(i, t)` and a small `ε > 0`.

The imputed value is:

```text
\hat{x}_{it} =
  [ Σ_{ℓ ∈ N_k(i,t)} α_{iℓt} x_{ℓt} ] /
  [ Σ_{ℓ ∈ N_k(i,t)} α_{iℓt} ]
```

### Categorical Imputation

For a categorical target variable `t`, let `v` range over the distinct categories observed among the selected donors. Define the weighted vote score:

```text
V_i(v) = Σ_{ℓ ∈ N_k(i,t)} α_{iℓt} 1[x_{ℓt} = v]
```

Then the imputed category is:

```text
\hat{x}_{it} = argmax_v V_i(v)
```

with deterministic tie-breaking.

### Single Imputation Operator

For one pass over all missing entries, define the imputation operator `T` such that:

```text
X^{(1)} = T(X)
```

where each missing entry is replaced by the corresponding numeric or categorical estimator computed from the observed donor set.

If iterative refinement is enabled in later versions, then:

```text
X^{(m+1)} = T(X^{(m)})
```

for iteration index `m`, but the default protocol is a single-pass engine.

### Multiple Imputation Operator

Let `M` be the number of completed datasets requested.

For each completed dataset `r ∈ {1, …, M}`, define a stochastic operator `T_r` that applies the same donor search but samples from the `k` nearest donors using normalized weights:

```text
P_{iℓt} = α_{iℓt} / [ Σ_{q ∈ N_k(i,t)} α_{iqt} ]
```

Then the `r`-th completed dataset is:

```text
X^{(r)} = T_r(X)
```

The distance model and neighbor set size are identical across imputations; randomness enters only through donor sampling within `N_k(i, t)`. The stochastic draw is repeated independently across `r`, subject to the same random seed policy. For numeric targets, the sampled donor value is copied directly. For categorical targets, the sampled donor category is copied directly.

The multiple-imputation output is the collection:

```text
{ X^{(1)}, X^{(2)}, …, X^{(M)} }
```

with associated provenance metadata.

### Algorithm 1: FM-KNN Imputation

```text
Input:
  - data matrix or data frame X
  - number of neighbors k
  - number of completed datasets M
  - variable weights w
  - scaling rule for numeric variables

Output:
  - one completed dataset if M = 1
  - a collection of M completed datasets if M > 1

Procedure:
  1. Detect variable types and construct the missingness mask R.
  2. Scale numeric variables using observed values only.
  3. For each target variable t with missing entries:
       a. Define the donor set D_t.
       b. For each receiver row i with R_it = 0:
            i.   Compute masked distances to all donors ℓ ∈ D_t.
            ii.  Select the k nearest finite-distance donors.
            iii. Compute the numeric or categorical estimator.
            iv.  Fill x_it in the single-imputation output.
  4. If M > 1:
       a. Repeat the donor selection step M times under the same distance model.
       b. Generate M completed datasets by stochastic donor sampling from the k nearest neighbors.
  5. Restore original variable types and return the completed object(s).
```

## Distance Rule

For numeric variables, the masked distance is:

```text
d(i, j) = sqrt(
  sum_p mask_ijp * w_p * (x_ip - x_jp)^2
  / sum_p mask_ijp * w_p
)
```

Where:

- `mask_ijp = 1` if variable `p` is observed in both rows, else `0`
- `w_p` is the variable weight
- the denominator prevents pairs with fewer shared variables from being unfairly favored or rejected

## Imputation Rule

### Numeric target

Use a distance-weighted mean of the `k` nearest donors:

```text
weight_j = 1 / (distance_j + epsilon)
```

The imputed value is the weighted average of donor target values.

### Categorical target

Use a distance-weighted majority vote among the `k` nearest donors.

If a tie occurs, resolve it deterministically.

## Version Scope

### v1.0

First release with the fast engine in place.

Supported scope:

- numeric data
- mixed data: numeric, factor, character, logical
- masked distance using jointly observed variables
- distance-weighted mean for numeric variables
- distance-weighted majority vote for categorical variables
- single and multiple imputation
- C++ backend via `cpp11`
- R wrapper API
- `complete()`
- `print()`
- `summary()`
- benchmark against `VIM::kNN()`
- optional imputation indicators

### v1.1

Follow-up performance and usability improvements:

- parallelization
- approximate nearest neighbors
- more distance options
- better diagnostics

## Planned Feature Roadmap

### v1.0

- C++ implementation
- performance-first optimization
- numeric + mixed data support
- partial top-`k` search
- weighted aggregation
- clean R API

### v1.1

- parallelization
- approximate search options if needed
- richer diagnostics

### v1.2

- integration with `mimar`
- simulation benchmarking
- publication-ready vignette

## Suggested User API

Minimal:

```r
imp <- missknn(
  data,
  k = 5,
  weights = "distance",
  scale = TRUE,
  max_iter = 1
)

completed <- complete(imp)
```

Extended:

```r
imp <- missknn(
  data,
  k = 7,
  distance = "masked",
  numeric = "weighted_mean",
  categorical = "weighted_mode",
  add_indicator = TRUE
)
```

## Internal Pipeline

```text
missknn()
  ├─ validate_data()
  ├─ detect_types()
  ├─ encode_data()
  ├─ scale_numeric()
  ├─ build_metadata()
  ├─ cpp_missknn_numeric()
  ├─ cpp_missknn_mixed()
  ├─ restore_types()
  └─ new_missknn_object()
```

C++ layer:

```text
src/
  cpp_missknn_numeric.cpp
  cpp_missknn_mixed.cpp
  masked_distance.cpp
  topk.cpp
  aggregate.cpp
```

## Implementation Principles

- Keep the default path simple and fast.
- Use `data.table` as the primary data structure for tabular handling in R.
- Separate numeric and mixed-data logic.
- Avoid recursive imputation in the core engine.
- Minimize data copies.
- Prefer matrix operations over row-by-row loops where possible.
- Keep dependencies minimal outside the C++ backend.
- Put distance computation, top-`k` neighbor search, and aggregation in C++ via `cpp11` in v1.

## Output Contract

The main result should be either a single completed dataset or a collection of multiply imputed datasets plus metadata sufficient for inspection and reproducibility.

The returned object should preserve:

- original data
- completed data
- multiple completed datasets when requested
- missingness indicators
- `k`
- distance settings
- scaling settings
- variable type information
- `data.table` compatibility for input and output workflows

## Benchmark Claim

The defensible claim is:

> `missknn` is a faster and simpler KNN imputation engine for tabular data.

Avoid claiming superiority in imputation quality over `VIM::kNN()` unless supported by empirical evaluation.

## Validation Plan

The package should be validated on:

- fully numeric datasets with artificial missingness
- mixed datasets after categorical support lands
- runtime comparison against `VIM::kNN()`
- memory usage comparison
- deterministic reproducibility checks

## Updated Abstract

`missknn` is a C++-backed R package for fast k-nearest neighbor imputation of missing values in tabular data. Unlike general-purpose imputation frameworks, `missknn` focuses on a compact and performance-oriented implementation of KNN imputation. The core algorithm uses masked distances computed only on jointly observed variables, partial top-`k` neighbor selection, distance-weighted aggregation, and a minimal user interface. The package is designed as a lightweight alternative for users who need fast, transparent, and reproducible single imputation and multiple imputation without the overhead of broader imputation frameworks.

## Documentation Plan

Required documentation artifacts:

- `README.md` with the main positioning statement
- this protocol document
- function reference for `missknn()`
- function reference for `complete()`
- benchmark vignette or report

## Acceptance Criteria for v1.0

Version 1.0 is acceptable when:

- the API is stable enough for early users
- numeric and mixed masked KNN work end-to-end
- single and multiple imputation work end-to-end
- output is reproducible
- documentation matches implementation
- the package can be benchmarked against `VIM::kNN()`

## Summary

`missknn` should be positioned as a modern, compact, performance-oriented KNN imputation engine.

Its main novelty is speed by design: the expensive distance, neighbor selection, and aggregation work lives in compiled C++, while R and `data.table` handle orchestration and output assembly.

`VIM::kNN()` remains the flexible reference implementation.

`missknn` wins by being:

- simpler
- smaller
- faster in the common path
- easier to explain
- easier to maintain
