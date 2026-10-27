# vgi4bio 0.0.0.9002

* New `fun_rel_abund()`: extracts the relative abundance of each species on
  each site from a fitted model, with delta-method confidence intervals,
  either rescaled to sum to one per species (`scale = "sum"`) or per unit of
  effort of the standardized dataset (`scale = "none"`).
* `fun_abund_estim()` gains a `standardized` argument naming the standardized
  dataset (known effort). Previously, the dataset fixed as reference was the
  first one in alphabetical order, which gave wrong estimates when the
  opportunistic dataset came first. Without `standardized`, a warning is now
  issued. The result also contains `standardized`.
* New long test reproducing Table 2 of Giraud et al. (2016)
  (run with `VGI4BIO_LONG_TESTS=true`).
* Vignette: extraction of relative abundances and comparison with STOC.

# vgi4bio 0.0.0.9001

* Package passes `R CMD check`: valid DESCRIPTION, complete documentation,
  documented `Aquitaine` data set, working tests.
* `fun_pearsons()` and `fun_braycurt()` are exported and take `newdata`
  explicitly.
* Warning for species with only zero counts in a dataset (structural zeros).
