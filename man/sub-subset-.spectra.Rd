% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.spectra.R
\name{[[<-.spectra}
\alias{[[<-.spectra}
\title{Assign values to spectra by wavelength}
\usage{
\method{[[}{spectra}(spectra, wavelength, j, values) <- value
}
\arguments{
\item{spectra}{Object of class \code{spectra}}

\item{wavelength}{Wavelength vector to select}

\item{j}{
    indices specifying elements to extract or replace.  Indices are
    \code{numeric} or \code{character} vectors or empty (missing) or
    \code{NULL}.  Numeric values are coerced to integer as by
    \code{\link[base]{as.integer}} (and hence truncated towards zero).
    Character vectors will be matched to the \code{\link[base]{names}} of the
    object (or for matrices/arrays, the \code{\link[base]{dimnames}}):
    see \sQuote{Character indices} below for further details.

    For \code{[}-indexing only: \code{i}, \code{j}, \code{\dots} can be
    logical vectors, indicating elements/slices to select.  Such vectors
    are recycled if necessary to match the corresponding extent.
    \code{i}, \code{j}, \code{\dots} can also be negative integers,
    indicating elements/slices to leave out of the selection.

    When indexing arrays by \code{[} a single argument \code{i} can be a
    matrix with as many columns as there are dimensions of \code{x}; the
    result is then a vector with elements corresponding to the sets of
    indices in each row of \code{i}.

    An index value of \code{NULL} is treated as if it were \code{integer(0)}.
  }

\item{values}{Vector or matrix of values to assign}
}
\description{
Assign values to spectra by wavelength
}
