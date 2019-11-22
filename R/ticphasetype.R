#' @details
#'
#' \code{ticphasetype} contains several functions for representing statistics that are
#' commonly used in population genomics. \code{\link{t_mrca}} and \code{\link{t_total}}
#' can be used for modelling the time until the most recent common ancestor and the
#' total tree length for Kingman's coalescence using a user-friendly interface. These
#' two quantities are represented using a continuous phase-type distribution
#' (the \code{\link{cont_phase_type}} class). The user can use the generator function
#' for this class if they want phase-type representations based on other coalescent
#' models.
#'
#' On the other hand, \code{ticphasetype} can also model the mutational process. The package
#' contains a generator function for the block counting process of Kingsman's coalescent
#' (\code{\link{kingsman}}). By reward-transforming this continuous phase-type representation,
#' \code{ticphasetype} creates discrete phase-type distributions of singletons, doubletons and
#' related statistics (\code{\link{itons}}), for the total number of segregating sites
#' (\code{\link{segsites}}), and for the tail statistic (\code{\link{tailstat}}). All these
#' quantities are represented using the \code{\link{disc_phase_type}} class, whose generator
#' function can also be used for user-tailored discrete phase-type distributions.
#'
#' The mean and the variance of the any phase-type representation can easily be computed
#' using \code{mean} and \code{\link{var}} respectively. Moreover, a readable
#' summary of a distribution can be obtained using the generic function \code{summary}.
#'
#' This package also contains functions for the density (\code{\link{dphtype}}), quantiles
#' (\code{\link{qphtype}}), distribution \code{\link{pphtype}} and random draw generator
#' (\code{\link{rphtype}}) for both discrete and continuous phase-type distributions.
#' Additionally, the user can calculate and plot the site-frequency spectrum using
#' \code{\link{sfs}}.
#'
#' For more details see the vignette by running \code{vignette("ticphasetype")}.
#'
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"
