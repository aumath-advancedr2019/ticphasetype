#' @details
#' This package contains two generator functions, namely \code{\link{cont_cont_phase_type}} and
#' \code{\link{disc_phase_type}} for continuous and discrete phase-type distributions
#' respectively. \code{cont_phase_type()} can easily be used for modelling the time until
#' the most recent common ancestor and the total tree length using a user-friendly
#' interface. On the other hand, \code{disc_phase_type()} generates phase-type
#' representations of i-tons (singletons, doubletons, etc.) and related statistics.
#' These functions also allow the user to build phase-type distributions at a higher
#' level by specifying their own sub-intensity matrix and initial probabilities.
#'
#' The mean and the variance of the phase-type distributions can easily be computed
#' using \code{mean} and \code{\link{var}} respectively. Moreover, a readable
#' summary of a distribution can be obtained using the generic function \code{summary}.
#'
#' This package also contains functions for the density (\code{\link{dphtype}}), quantiles
#' (\code{\link{qphtype}}), distribution \code{\link{pphtype}} and random draw generator
#' (\code{\link{rphtype}}) for the phase-type distributions. Additionally, the user can
#' calculate and plot the site-frequency spectrum using \code{\link{sfs}}.
#'
#' For more details see the vignette by running \code{vignette("ticphasetype")}.
#'
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"

#> [1] "_PACKAGE"
