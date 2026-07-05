# missknn Benchmark Report

This benchmark uses 30 simulations, MCAR amputation from `mimar`, linear regression coefficients, and `bench::mark` for runtime.

## Figures

![Runtime plot](runtime_plot.png)

![Bias plot](bias_plot.png)

![MSE plot](mse_plot.png)

## Runtime Table

Table: Mean runtime across 30 simulations

|method     | runtime_sec| sd_runtime|
|:----------|-----------:|----------:|
|missknn    |      0.0043|     0.0011|
|missForest |      0.5130|     0.1757|

## Method-Level Error Table

Table: Mean coefficient error across 30 simulations

|method     | abs_bias|    mse|
|:----------|--------:|------:|
|missknn    |   0.0335| 0.0018|
|missForest |   0.0394| 0.0024|

## Coefficient-Level Table

Table: Coefficient-level absolute bias and MSE

|method     |term        | abs_bias|    mse|
|:----------|:-----------|--------:|------:|
|missknn    |(Intercept) |   0.0270| 0.0011|
|missForest |(Intercept) |   0.0278| 0.0011|
|missknn    |x1          |   0.0303| 0.0014|
|missForest |x1          |   0.0428| 0.0027|
|missknn    |x2          |   0.0347| 0.0019|
|missForest |x2          |   0.0387| 0.0020|
|missknn    |x3          |   0.0347| 0.0020|
|missForest |x3          |   0.0410| 0.0025|
|missknn    |x4          |   0.0329| 0.0017|
|missForest |x4          |   0.0409| 0.0029|
|missknn    |x5          |   0.0412| 0.0026|
|missForest |x5          |   0.0450| 0.0029|

