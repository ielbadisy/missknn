# missknn Real-Data Benchmark Report

Metrics follow Stekhoven & Buhlmann (2012): NRMSE for numeric variables
and PFC (proportion falsely classified) for categorical variables, both
pooled across every artificially-masked entry of that type rather than
averaged per column first. 20% MCAR imposed on complete-case ground
truth for each dataset.

## Dataset selection

Four biomedical/clinical datasets from `biostatlab`, kept to a size
appropriate for an Application Note while spanning a range of outcomes:
`heart_failure` (missknn wins outright), `crc_mondaca2020` (missknn wins
on speed and numeric accuracy, loses on categorical accuracy),
`arthritis` (missknn wins on speed only), and `metabric_clinical`
(missForest fails outright on a rare categorical level; missknn and
missRanger both complete, missRanger slightly more accurate).
`diabetes_prediction` (n = 100,000) is reported separately as a real-data
scaling point.

## Figures

![Runtime by dataset](real_data_speed_plot.png)

![NRMSE by dataset](real_data_nrmse_plot.png)

## Main Results Table

Table: Runtime, NRMSE, and PFC on real biomedical datasets (20% MCAR)

|dataset           |method     |    n|  p| runtime_sec|  nrmse|    pfc|
|:-----------------|:----------|----:|--:|-----------:|------:|------:|
|heart_failure     |missknn    |  299| 13|       0.050| 0.3546|     NA|
|heart_failure     |missForest |  299| 13|       2.193| 0.3769|     NA|
|heart_failure     |missRanger |  299| 13|       0.068| 0.3867|     NA|
|heart_failure     |VIM::kNN   |  299| 13|       0.263| 0.3892|     NA|
|crc_mondaca2020   |missknn    |  457| 11|       0.057| 0.4081| 0.4346|
|crc_mondaca2020   |missForest |  457| 11|       4.457| 0.4220| 0.4191|
|crc_mondaca2020   |missRanger |  457| 11|       0.777| 0.4199| 0.3969|
|crc_mondaca2020   |VIM::kNN   |  457| 11|       0.512| 0.5651| 0.4435|
|arthritis         |missknn    | 4856| 11|       0.468| 0.7928| 0.2302|
|arthritis         |missForest | 4856| 11|      28.617| 0.7115| 0.2661|
|arthritis         |missRanger | 4856| 11|       5.006| 0.7164| 0.2197|
|arthritis         |VIM::kNN   | 4856| 11|      26.443| 0.8365| 0.2459|
|metabric_clinical |missknn    | 1839| 11|       0.451| 0.5821| 0.2983|
|metabric_clinical |missForest | 1839| 11|          NA|     NA|     NA|
|metabric_clinical |missRanger | 1839| 11|       3.497| 0.5725| 0.2898|
|metabric_clinical |VIM::kNN   | 1839| 11|       3.475| 0.6473| 0.3362|

## Real-Data Scaling Point (n = 100,000)

Table: diabetes_prediction (n = 100,000): missForest excluded as computationally intractable at this n

|dataset             |method     |     n|  p| runtime_sec|  nrmse| pfc|
|:-------------------|:----------|-----:|--:|-----------:|------:|---:|
|diabetes_prediction |missknn    | 1e+05|  9|      29.931| 0.3471|  NA|
|diabetes_prediction |missRanger | 1e+05|  9|      73.406| 0.3616|  NA|

