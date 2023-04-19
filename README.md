# ecurve

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

[![codecov](https://codecov.io/gh/phac-nml-phrsd/ecurve/branch/main/graph/badge.svg?token=HKCD5AE7KP)](https://codecov.io/gh/phac-nml-phrsd/ecurve)

Enhanced Standard Curve Model

This R library implements a statistical model ([Schmidt et al. 2022](https://www.frontiersin.org/articles/10.3389/fmicb.2023.1048661/full) )  that provides an improved representation of qPCR assays standard curve data at low concentrations while converging asymptotically upon conventional log-linear regression at high concentrations.

See the introductory vignette shipped with the package for more details.

## Installation

Use the `devtools` library to install `ecurve`: 

`devtools::install_github(repo="phac-nml-phrsd/ecurve", build_vignettes = TRUE)`

