README
================
Karine Princé
12/11/2020

<!-- README.md is generated from README.Rmd. Please edit that file -->

# vgi4bio

<!-- badges: start -->
<!-- badges: end -->

The package vgi4bio contains functions to estimate relative species
abundances from opportunistic records, using the data combination
methods of Giraud et al. (2016, *Biometrics*,
[doi:10.1111/biom.12431](https://doi.org/10.1111/biom.12431)). The
estimation of habitat preferences (Coron et al. 2018, *Environmental and
Ecological Statistics*,
[doi:10.1007/s10651-018-0398-2](https://doi.org/10.1007/s10651-018-0398-2))
is under development. To have an overview of the package, read the
vignette after installing the package by:

``` r
vignette("tutorial", package = "vgi4bio")
```

## Installing *vgi4bio*

To install the development version from GitHub (with the vignette):

1.  Install the release version of devtools from CRAN with
    install.packages(“devtools”).

2.  Make sure you have a working development environment.

- Windows: Install Rtools.
- Mac: Install Xcode from the Mac App Store.
- Linux: Install a compiler and various development libraries (details
  vary across different flavors of Linux).

Then:

``` r
library(devtools)
install_github("KarinePrince/vgi4bio", build_vignettes = TRUE)
```

The package is not (yet) available on CRAN.

Once installed, the package can be loaded using:

``` r
library(vgi4bio)
```

*The package requires R \>= 3.5.*
