# missknn Real-Data Benchmark Report

Five real biomedical/clinical datasets (from the `biostatlab` package) with
20% MCAR imposed column-wise on complete-case ground truth, comparing
`missknn`, `missMDA::imputePCA` (numeric columns only), `missForest`, and
`VIM::kNN` (numeric columns only):

- **pbc** (n = 276): Mayo Clinic primary biliary cirrhosis trial data, all
  numeric/integer-coded.
- **heart_failure** (n = 299): UCI heart failure clinical records, all
  numeric/integer-coded.
- **framingham** (n = 5209): larger real cardiovascular cohort, all
  numeric/integer-coded.
- **metabric_clinical** (n complete cases, p = 11): breast cancer clinical
  variables from the METABRIC cohort, mixing 6 numeric and 5 categorical
  (factor) variables -- `missMDA::imputePCA` and `VIM::kNN` are only
  applied to the numeric subset here since they do not support factors
  directly (`missMDA` needs the separate `imputeFAMD` routine for mixed
  data); `missknn` and `missForest` impute numeric and categorical targets
  from the same call.
- **crc_mondaca2020** (n complete cases, p = 11): metastatic colorectal
  cancer clinical/genomics cohort (Mondaca et al. 2020), mixing 6 numeric
  (including tumor mutational burden and MSI score) and 5 categorical
  variables.

## Results Table

Table: Runtime, numeric MSE, and categorical accuracy on real biomedical datasets (20% MCAR)

|dataset           |method     | runtime_sec|  numeric_mse| factor_accuracy|
|:-----------------|:----------|-----------:|------------:|---------------:|
|pbc               |missknn    |       0.095| 3.143302e+05|              NA|
|pbc               |missMDA    |       0.087| 1.958058e+05|              NA|
|pbc               |missForest |       2.756| 2.183684e+05|              NA|
|pbc               |VIM::kNN   |       0.354| 2.522463e+05|              NA|
|heart_failure     |missknn    |       0.028| 1.002789e+08|              NA|
|heart_failure     |missMDA    |       0.028| 1.193960e+08|              NA|
|heart_failure     |missForest |       1.480| 1.190990e+08|              NA|
|heart_failure     |VIM::kNN   |       0.326| 1.280750e+08|              NA|
|framingham        |missknn    |       5.562| 2.419981e+02|              NA|
|framingham        |missMDA    |       0.511| 2.644414e+02|              NA|
|framingham        |missForest |     212.981| 2.394146e+02|              NA|
|framingham        |VIM::kNN   |      59.431| 3.811433e+02|              NA|
|metabric_clinical |missknn    |       0.632| 2.068043e+02|          0.9380|
|metabric_clinical |missMDA    |       0.060| 1.921224e+02|              NA|
|metabric_clinical |missForest |          NA|           NA|              NA|
|metabric_clinical |VIM::kNN   |       1.603| 2.919995e+02|              NA|
|crc_mondaca2020   |missknn    |       0.050| 1.155001e+01|          0.9112|
|crc_mondaca2020   |missMDA    |       0.361| 1.296953e+01|              NA|
|crc_mondaca2020   |missForest |       4.501| 1.334407e+01|          0.9068|
|crc_mondaca2020   |VIM::kNN   |       0.159| 1.829211e+01|              NA|

## Note on the missForest failure

`missForest` failed outright on `metabric_clinical` (`NA` above): one of its
categorical predictors has a rare level, and if 20% MCAR happens to mask
every instance of that level in a given draw, `missForest`'s internal
per-level out-of-bag error accounting has no observations left to compare
against and errors out. `missknn`'s masked distance has no such failure
mode -- a rare level is simply a rare donor category, not a special case.

