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
|missknn    |      0.0044|     0.0012|
|missForest |      0.5333|     0.1842|
|missMDA    |      0.0076|     0.0036|
|VIM::kNN   |      0.0611|     0.0121|

## Method-Level Error Table

Table: Mean coefficient error across 30 simulations

|method     | abs_bias|    mse|
|:----------|--------:|------:|
|missknn    |   0.0335| 0.0018|
|missForest |   0.0394| 0.0024|
|missMDA    |   0.0358| 0.0020|
|VIM::kNN   |   0.0390| 0.0024|

## Coefficient-Level Table

Table: Coefficient-level absolute bias and MSE

|method     |term        | abs_bias|    mse|
|:----------|:-----------|--------:|------:|
|missknn    |(Intercept) |   0.0270| 0.0011|
|missForest |(Intercept) |   0.0278| 0.0011|
|missMDA    |(Intercept) |   0.0269| 0.0011|
|VIM::kNN   |(Intercept) |   0.0307| 0.0012|
|missknn    |x1          |   0.0303| 0.0014|
|missForest |x1          |   0.0428| 0.0027|
|missMDA    |x1          |   0.0339| 0.0016|
|VIM::kNN   |x1          |   0.0407| 0.0026|
|missknn    |x2          |   0.0347| 0.0019|
|missForest |x2          |   0.0387| 0.0020|
|missMDA    |x2          |   0.0354| 0.0019|
|VIM::kNN   |x2          |   0.0463| 0.0033|
|missknn    |x3          |   0.0347| 0.0020|
|missForest |x3          |   0.0410| 0.0025|
|missMDA    |x3          |   0.0373| 0.0024|
|VIM::kNN   |x3          |   0.0445| 0.0028|
|missknn    |x4          |   0.0329| 0.0017|
|missForest |x4          |   0.0409| 0.0029|
|missMDA    |x4          |   0.0389| 0.0024|
|VIM::kNN   |x4          |   0.0380| 0.0024|
|missknn    |x5          |   0.0412| 0.0026|
|missForest |x5          |   0.0450| 0.0029|
|missMDA    |x5          |   0.0425| 0.0028|
|VIM::kNN   |x5          |   0.0341| 0.0020|

