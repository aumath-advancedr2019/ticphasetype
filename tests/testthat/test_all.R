

#' @import testthat


# 1 Test continuous phase type functions

test_that("dphtype gives sensible results", {
  pht = t_mrca(6)
  observed_1 = dphtype(seq(0, 10), pht)
  expected_1 = c(
    0,
    0.5338932351099173745368,
    0.276755865491356389807,
    0.1060255272732280595882,
    0.0392148823841512172983,
    0.01443681909321984507821,
    0.005311530217966157346732,
    0.001954028721245163000703,
    0.0007188482861231658633053,
    0.0002644495701168376252006,
    9.72855632754532318688e-05
  )
  expect_equal(observed_1, expected_1)
})


test_that("qphtype gives sensible results", {
  pht = t_mrca(6)
  observed_1 = qphtype(seq(0, 0.999, 0.1),
                       pht)
  expected_1 = c(
    0.0000000000000000000000,
    0.6072698798036705314374,
    0.8147133568055460184354,
    1.0027356505258440133588,
    1.1945201743795064164289,
    1.4040717057458642624823,
    1.6471119434939456294131,
    1.9490937401898580372261,
    2.3641840256832784561425,
    3.0628993748782011863341
  )
  expect_equal(observed_1, expected_1, tolerance = 1e-05)

})


test_that("pphtype gives sensible results", {
  pht = t_mrca(6)
  observed_1 = pphtype(seq(0, 1, 0.1), pht)
  expected_1 = c(
    0,
    0.0001280113094311863264352,
    0.002417571282576025382127,
    0.01120471786207932751722,
    0.02971417578333634956778,
    0.05867961149384781638361,
    0.09691010576874470316966,
    0.1422592616731122028284,
    0.1923621809256920167641,
    0.2450517395979424639663,
    0.2985364854629098951833
  )
  expect_equal(observed_1, expected_1)

})



test_that("rphtype gives sensible results", {
  pht = t_mrca(6)
  rph = rphtype(10000, pht)
  expect_equal(mean(rph), 1.66, tolerance = 0.2)
  expect_equal(stats::var(rph), 1.14, tolerance = 0.5)

})


test_that("reward transformations is consistent", {
  n = 5

  obj = t_mrca(n)
  observed = rewardtransformation((1):(n-1), obj$init_probs, obj$subint_mat)$subint_mat

  expected = matrix(c(-10,10,0,0.00,
                      0,-3,3,0.00,
                      0,0,-1,1.00,
                      0,0,0,-0.25), ncol = n-1, nrow = n-1, byrow = T)


  expect_equal(expected, observed, tolerance = 1e-6)


  n = 8

  obj = t_mrca(n)
  observed = rewardtransformation((1):(n-1), obj$init_probs, obj$subint_mat)$subint_mat

  expected = matrix(c(-28, 28.0, 0.0, 0.0, 0.0, 0.0, 0.0000000,
                      0, -10.5, 10.5, 0.0, 0.0, 0.0, 0.0000000,
                      0, 0.0, -5.0, 5.0, 0.0, 0.0, 0.0000000,
                      0, 0.0, 0.0, -2.5, 2.5, 0.0, 0.0000000,
                      0, 0.0, 0.0, 0.0, -1.2, 1.2, 0.0000000,
                      0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5000000,
                      0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1428571), ncol = n-1, nrow = n-1, byrow = T)


  expect_equal(expected, observed, tolerance = 1e-6)
})


test_that("RewTransform wrapper for rewardtransformation is consistent", {
  n = 5

  mphobj = kingsman(n)

  reward_length = length(mphobj$init_probs)
  transformed = RewTransform(mphobj, 0:(reward_length-1), 2)
  transformed$subint_mat

  expected = matrix(c(0.1428571, 0.1714286, 0.2142857, 0.1942857, 0.1904762,
                      0.0000000, 0.4000000, 0.0000000, 0.3200000, 0.1666667,
                      0.0000000, 0.0000000, 0.5000000, 0.1333333, 0.2777778,
                      0.0000000, 0.0000000, 0.0000000, 0.8000000, 0.0000000,
                      0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.8333333), nrow = 5, ncol = 5, byrow = TRUE)

  expect_equal(expected, transformed$subint_mat, tolerance = 1e-6)

})


# clean data
#
# library(tidyverse)
# data1 = read_delim(
#   "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf",
#   "\t",
#   escape_double = FALSE,
#   comment = "##",
#   trim_ws = TRUE
# ) %>% select(starts_with(c("HG", "NA")))
# #data1 = data1[10:ncol(data1), 1:1000]
# data1 %>% select(
#   yoruban = NA18881,
#   iberian = HG01620,
#   chinese = HG00513,
#   peruvian = HG02090,
#   luhya = NA19019
# )

# compute the SFS
