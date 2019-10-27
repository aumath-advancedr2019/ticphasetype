


# This block makes a plot that shows that the dphtype and rphtype agree.
sanity_check_rphdph = function(OTUs = 5, draws = 1000000) {




  # rphtype data
  r_data = tibble(x = rphtype(draws, OTUs),
                  p = rep(NA, length(x)),
                  data = "random sampling")# %>% ggplot(aes(x)) + geom_histogram(binwidth = 0.01)


  # dphtype data
  d_data = tibble(x = seq(0, 10, 0.01),
                  p = dphtype(seq(0, 10, 0.01), generate_init_row(OTUs-1), generate_subint_mat(OTUs)),
                  data = "theoretical")
  #
  ggplot(mapping = aes(color = data)) +
    geom_point(mapping = aes(x, p), data = d_data) +
    geom_line(mapping = aes(r_data$x), data = r_data, stat = "density", size = 1)
}


sanity_check_rphdph(draws = 1000000)
