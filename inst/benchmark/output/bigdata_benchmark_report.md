# missknn Big-Data Benchmark Report

Runtime and imputation accuracy vs. sample size, from n = 1,000 up to
n = 100,000 (5 numeric columns, 20% MCAR per column). Each (n, method) run
executes in its own fresh `Rscript` process, so timings are not distorted
by heap growth accumulated across earlier, larger runs in a shared session.
`missForest` is only run up to n = 5,000 since its per-tree cost makes
larger n impractical; `missknn` and `missMDA` run across the full grid.

## Takeaway

`missknn` and `missMDA` are essentially tied on accuracy (MSE) at every n;
both dramatically outperform `missForest`, whose runtime already becomes
impractical at n = 5,000. On speed, `missMDA`'s lower fixed overhead keeps
it faster at small/medium n, but `missknn`'s cost grows more slowly with n
(capped-donor masked-KNN search is close to linear, vs. `missMDA`'s
iterative PCA passes), so `missknn` overtakes `missMDA` and pulls ahead as
n grows into genuine big-data territory (roughly n >= 100,000 here).

## Figures

![Runtime plot](bigdata_runtime_plot.png)

![MSE plot](bigdata_mse_plot.png)

## Results Table

Table: Runtime and MSE by method and sample size (one fresh Rscript process per run)

|     n|method     | runtime_sec|     mse|
|-----:|:----------|-----------:|-------:|
| 1e+03|missknn    |      0.0496| 0.18962|
| 1e+03|missMDA    |      0.0101| 0.19400|
| 1e+03|missForest |      4.3177| 0.21800|
| 5e+03|missknn    |      0.1777| 0.19814|
| 5e+03|missMDA    |      0.0351| 0.19823|
| 5e+03|missForest |     73.9245| 0.21293|
| 2e+04|missknn    |      0.1836| 0.19876|
| 2e+04|missMDA    |      0.1673| 0.19903|
| 1e+05|missknn    |      0.2503| 0.19955|
| 1e+05|missMDA    |      0.9757| 0.19955|

