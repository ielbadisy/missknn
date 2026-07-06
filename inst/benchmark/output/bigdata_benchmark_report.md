# missknn Big-Data Benchmark Report

Runtime and imputation accuracy (NRMSE) vs. sample size, from n = 1,000
up to n = 100,000 (5 numeric columns, 20% MCAR per column). Each
(n, method) run executes in its own fresh `Rscript` process, so timings
are not distorted by heap growth accumulated across earlier, larger runs
in a shared session. `missForest` and `VIM::kNN` are only run up to
n = 5,000 (per-tree cost and O(n^2) distance search respectively make
larger n impractical); `missknn` and `missRanger` run across the full
grid.

## Takeaway

`missForest`'s runtime already becomes impractical at n = 5,000 (over a
minute) and would take hours at n = 100,000; `VIM::kNN`'s distance search
is similarly infeasible past a few thousand rows. `missknn`'s
capped-donor masked-KNN search stays close to linear in n and keeps
running in well under a second per column set even at n = 100,000, which
is the practical point of the donor-cap design: search cost does not
grow with the full donor pool once n exceeds the cap. `missRanger` also
scales to n = 100,000 but at substantially higher runtime.

## Figures

![Runtime plot](bigdata_runtime_plot.png)

![NRMSE plot](bigdata_mse_plot.png)

## Results Table

Table: Runtime and NRMSE by method and sample size (one fresh Rscript process per run)

|     n|method     | runtime_sec|   nrmse|
|-----:|:----------|-----------:|-------:|
| 1e+03|missknn    |      0.0476| 1.00223|
| 1e+03|missRanger |      0.0942| 1.15781|
| 1e+03|missForest |      3.1344| 1.07461|
| 1e+03|VIM::kNN   |      0.4911| 1.16132|
| 5e+03|missknn    |      0.1337| 1.00020|
| 5e+03|missRanger |      1.0260| 1.14478|
| 5e+03|missForest |     73.0740| 1.03685|
| 5e+03|VIM::kNN   |     10.3536| 1.15992|
| 2e+04|missknn    |      0.1013| 0.99998|
| 2e+04|missRanger |      5.8439| 1.14941|
| 1e+05|missknn    |      0.2006| 1.00000|
| 1e+05|missRanger |     48.0944| 1.14435|

