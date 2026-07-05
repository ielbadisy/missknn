# missknn Correlated- and Mixed-Data Validation

## Study 1: correlated numeric data

n = 300, p = 6 numeric variables drawn from an equicorrelated Gaussian
(rho = 0.6), 20% MCAR per column.

|method     |     mse|
|:----------|-------:|
|missknn    | 0.12024|
|missForest | 0.11099|

## Study 2: mixed numeric + categorical data

n = 400, 3 correlated numeric variables plus 1 binary factor driven by one
of them, 20% MCAR per column. `missknn` and `missForest` handle numeric and
categorical targets from the same call.

|method     | numeric_mse| factor_accuracy|
|:----------|-----------:|---------------:|
|missknn    |     0.15832|          0.9700|
|missForest |     0.15441|          0.9775|

