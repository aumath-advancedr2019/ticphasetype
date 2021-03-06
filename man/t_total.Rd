% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/t_total.R
\name{t_total}
\alias{t_total}
\title{Total tree length}
\usage{
t_total(n)
}
\arguments{
\item{n}{positive integer larger than 1.}
}
\description{
This function creates a continuous phase-type representation of the total tree length
for a specified number of sequences and following Kingsman's coalescent model.
}
\examples{
ph_total <- t_total(4)
summary(ph_total)

}
