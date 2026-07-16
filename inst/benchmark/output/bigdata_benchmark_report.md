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
|  1,000|missknn    |      0.0131| 1.00223|
|  1,000|missRanger |      0.1630| 1.15781|
|  1,000|missForest |      3.1537| 1.07461|
|  1,000|VIM::kNN   |      0.6088| 1.16132|
|  5,000|missknn    |      0.0321| 1.00020|
|  5,000|missRanger |      1.3632| 1.14478|
|  5,000|missForest |     75.1881| 1.03685|
|  5,000|VIM::kNN   |     10.9076| 1.15992|
| 20,000|missknn    |      0.0456| 0.99998|
| 20,000|missRanger |     10.3372| 1.14941|
|100,000|missknn    |      0.0938| 1.00000|
|100,000|missRanger |     85.6539| 1.14435|

