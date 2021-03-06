% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sfs.R
\name{sfs}
\alias{sfs}
\title{Site Frequency Spectrum}
\usage{
sfs(n, theta = 2, plot = FALSE)
}
\arguments{
\item{n}{sample size}

\item{theta}{population mutation parameter}

\item{plot}{logical}
}
\description{
This function returns a list of expected counts of itons up to n-1, or a plot if specified.
It builds the site frequency spectrum using Kingman's coalescent. It works for different
values of `theta` (population mutation parameter).
}
\details{
By default, a list of expected values for each i-ton and its variance is returned.
If \code{plot=TRUE}, then a plot is returned.
}
\examples{
sfs(n = 5, theta = 3)

}
