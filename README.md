README
================
Karine Princé
12/11/2020

<!-- README.md is generated from README.Rmd. Please edit that file -->

# vgi4bio

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/KarinePrince/vgi4bio.svg?branch=master)](https://travis-ci.com/KarinePrince/vgi4bio)
<!-- badges: end -->

The package of vgi4bio contains functions to estimate relative species
abundances from opportunistic records, using data combination methods.
To have an overview of the package, read the vignette after installing
the package by:

``` r
vignette("tutorial", package = "vgi4bio")
```

## Installing *vgi4bio*

To install the development version from GitHub:

1.  Install the release version of devtools from CRAN with
    install.packages(“devtools”).

2.  Make sure you have a working development environment.

<!-- end list -->

  - Windows: Install Rtools.
  - Mac: Install Xcode from the Mac App Store.
  - Linux: Install a compiler and various development libraries (details
    vary across different flavors of Linux).

Then:

``` r
library(devtools)
install_github("KarinePrince/vgi4bio")
```

The stable version can be installed from
[CRAN](https://CRAN.R-project.org) using:

``` r
install.packages("vgi4bio")
```

Once installed, the package can be loaded using:

``` r
library(vgi4bio)
```

*Note that this package was build with R 3.6, so you might have to
update your R installation.*
