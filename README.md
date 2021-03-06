
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PEcAn Radiative Transfer Modeling module – standalone version

**Corresponding author**:

Alexey Shiklomanov

NASA Goddard Space Flight Center

<alexey.shiklomanov@nasa.gov>

## Installation

Easiest way to install is via `install_github` from the `devtools`
package.

``` r
devtools::install_github("ashiklom/pecanrtm")
```

If you want a specific branch, do `install_github(..., ref="<branch>")`.

From there, you should be able to load the package in your typical R
session.

*NOTE: On some OS X systems, the automatic specification of `gfortran`
will fail, causing the installation to abort. To fix this issue, you
will need to add something like the following the `~/.R/Makevars` file.*

    FC = gfortran

For more information, see the vignette
(`vignettes/pecanrtm.vignettes.Rmd`).

## Usage

Simulate leaf reflectance and transmittance with different versions of
the PROSPECT model.

``` r
library(PEcAnRTM)
#> 
#> Attaching package: 'PEcAnRTM'
#> The following object is masked from 'package:graphics':
#> 
#>     matplot
p4 <- prospect(c(N = 1.4, Cab = 40, Cw = 0.01, Cm = 0.005), version = "4")
p5 <- prospect(c(N = 1.4, Cab = 40, Car = 8, Cw = 0.01, Cm = 0.005), version = "5")
p5b <- prospect(c(
  N = 1.4, Cab = 40, Car = 8,
  Cbrown = 0.1, Cw = 0.01, Cm = 0.005
), version = "5B")
pd <- prospect(c(
  N = 1.4, Cab = 40, Car = 8, Canth = 10,
  Cbrown = 0.1, Cw = 0.01, Cm = 0.005
), version = "D")

wl <- seq(400, 2500)

plot(0, 0, type = "n", xlim = c(400, 2500), ylim = c(0, 1),
     xlab = "Wavelength (nm)", ylab = "Reflectance or 1 - Transmittance")
lines(wl, p4[, 1], col = 2)
lines(wl, 1 - p4[, 2], col = 2, lty = 2)
lines(wl, p5[, 1], col = 3)
lines(wl, 1 - p5[, 2], col = 3, lty = 2)
lines(wl, p5b[, 1], col = 4)
lines(wl, 1 - p5b[, 2], col = 4, lty = 2)
lines(wl, pd[, 1], col = 5)
lines(wl, 1 - pd[, 2], col = 5, lty = 2)
legend("topright", c("4", "5", "5B", "D"), col = 2:5, lty = 1)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

Simulate canopy reflectance with the PRO4SAIL model (PROSPECT 5B coupled
to the 4SAIL canopy RTM).

``` r
canopy <- pro4sail(c(
  N = 1.4, Cab = 40, Car = 8, Cbrown = 0, Cw = 0.01, Cm = 0.009,
  # Leaf angle distribution function parameters
  LIDFa = -0.35, LIDFb = -0.15, TypeLIDF = 1,
  # Leaf area index
  LAI = 3,
  # Hot spot parameter
  q = 0.01,
  # Solar zenith, degrees
  tts = 30,
  # Observer zenith, degrees
  tto = 10,
  # Sun-observer azimuth
  psi = 0,
  # Soil moisture fraction (0 = dry, 1 = wet)
  psoil = 0.7
))

matplot(canopy)
legend("topright", c("BHR", "HDR", "DHR", "BDR"), col = 1:4, lty = 1:4)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
