# PEcAn Radiative Transfer Modeling module -- standalone version

**Corresponding author**

Alexey Shiklomanov

NASA Goddard Space Flight Center

alexey.shiklomanov@nasa.gov

## Installation

Easiest way to install is via `install_github` from the `devtools` package.

```R
devtools::install_github("ashiklom/pecanrtm")
```

If you want a specific branch, do `install_github(..., ref="<branch>")`.

From there, you should be able to load the package in your typical R session.

*NOTE: On some OS X systems, the automatic specification of `gfortran` will 
fail, causing the installation to abort. To fix this issue, you will need to 
add something like the following the `~/.R/Makevars` file.*

```
FC = gfortran
```

For more information, see the vignette (`vignettes/pecanrtm.vignettes.Rmd`).
