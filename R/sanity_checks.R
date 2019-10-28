#' We ought to have a test-framework. But until we make one, I thought we could put some sanity checks into this file.


#' Sanity checks, showing that everything works as intended
#'
#' @param test unknown
#' @export



#' @import tibble
#' @import ggplot2


#' @describeIn sanity_check_rphdph
#'
#' Makes a simple comparison between theoretical and randomly sampled data from dphtype() and rphtype, respectively
#'
#' @param OTUs number of operational taxonomical units
#' @param draws number of random draws from rphtype wanted
#' @param granularity the difference between the numbers put into dphtype()
#'
#' @usage sanity_check_rphdph = function(OTUs, draws, granularity)
#'
#' @export


# This block makes a plot that shows that the dphtype and rphtype agree.
sanity_check_rphdph = function(OTUs = 5, draws = 1000000, granularity = 0.01) {

  # rphtype data
  init_probs = generate_init_row(OTUs - 1)
  subint_mat = generate_subint_mat(OTUs)
  r_data = tibble(x = rphtype(draws, init_probs, subint_mat, granularity),
                  p = rep(NA, length(x)),
                  data = "random sampling")# %>% ggplot(aes(x)) + geom_histogram(binwidth = granularity)


  # dphtype data
  d_data = tibble(x = seq(0, 10, granularity),
                  p = dphtype(seq(0, 10, granularity), generate_init_row(OTUs-1), generate_subint_mat(OTUs)),
                  data = "theoretical")
  #
  ggplot(mapping = aes(color = data)) +
    geom_point(mapping = aes(x, p), data = d_data) +
    geom_line(mapping = aes(r_data$x), data = r_data, stat = "density", size = 1) +
    labs(subtitle = paste0("Comparison of dphtype() and rphtype() with ", draws, " random draws and a gran. of ", granularity))
}

#sanity_check_rphdph(draws = 100000)
