## Update submission (1.1.2 -> 1.2.0)

`missknn()` now returns the completed data directly (a data frame, or a list
of completed data frames when `m >= 2`). The `missknn` object class and
`complete()` are removed; analysis and pooling are left to dedicated packages
such as mimar. `data.table` is replaced by `basetable`, only `NA` cells are
imputed, and ordered factors are kept. An `inst/CITATION` file is added.
See NEWS.md.

## Reverse dependencies

mimar (Imports) is the only reverse dependency. mimar 1.4.2, now on CRAN,
handles both the old and the new return value; its checks pass with this
version.

## Test environments

* Local: Ubuntu 24.04, R 4.5.1

## R CMD check results

0 errors | 0 warnings | 1 note

* "Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'". This flag comes from the local Ubuntu R
  toolchain (`Makeconf`), not from the package's `src/Makevars`.
