% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disc_phase_type_separated.R
\name{RewTransform}
\alias{RewTransform}
\title{Discretize MPH object with a reward vector}
\usage{
RewTransform(mph_obj, rewards, theta)
}
\arguments{
\item{mph_obj}{Multivariate Phase Type object generated either from mult_phase_type or kingsman function}

\item{rewards}{vector of non negative numbers}

\item{theta}{mutation parameter (positive number)}
}
\value{
A `disc_phase_type` object containing subintensity matrix (P), vector of initial probabilities (alpha) and defect (probability of not entering any transient
state prior to absorption)
}
\description{
Reward transformation of continuous Phase Type distribution into Discrete Phase Type
}
\examples{
RewTransform(kingsman(4), c(4, 2, 1, 0), 2)

}