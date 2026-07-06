# missknn Correlated- and Mixed-Data Validation

Metrics: NRMSE (numeric) and PFC (categorical), as in
Stekhoven & Buhlmann (2012); see inst/benchmark/benchmark_real.R.

## Study 1: correlated numeric data

n = 300, p = 6 numeric variables drawn from an equicorrelated Gaussian
(rho = 0.6), 20% MCAR per column.

|method     |   nrmse|
|:----------|-------:|
|missknn    | 0.82242|
|missForest | 0.80700|
|missRanger | 0.81623|
|VIM::kNN   | 0.88845|

## Study 2: mixed numeric + categorical data

n = 400, 3 correlated numeric variables plus 1 binary factor driven by one
of them, 20% MCAR per column.

|method     |   nrmse|    pfc|
|:----------|-------:|------:|
|missknn    | 0.85508| 0.2000|
|missForest | 0.91550| 0.2250|
|missRanger | 0.82490| 0.2000|
|VIM::kNN   | 0.95601| 0.2625|

