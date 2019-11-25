[![Build Status](https://travis-ci.com/aumath-advancedr2019/ticphasetype.svg?branch=master)](https://travis-ci.com/aumath-advancedr2019/ticphasetype)

# ticphasetype <img src="man/figures/TIC_logo.svg" align="right" />

This package implements core functions from phase-type theory that are useful in the analysis of biological data. It includes easy-to-use functions for modelling common statistics in evolutionary genetics, as well as a flexible high-level interface for more advanced users. You can learn more about the package in `vignette('ticphasetype')`.

## Details

ticphasetype contains several functions for representing statistics that are commonly used in population genomics. `t_mrca()` and `t_total()` can be used for modelling the time until the most recent common ancestor and the total tree length for Kingman's coalescence using a user-friendly interface. These two quantities are represented using a discrete phase-type distribution (the `cont_phase_type` class). The user can use the generator function for this class if they want phase-type representations based on other coalescent models.

On the other hand, ticphasetype can also model the mutational process. The package contains a generator function for the block counting process of Kingsman's coalescent (`kingsman()`). By reward-transforming this continuous phase-type representation, ticphasetype creates discrete phase-type distributions of singletons, doubletons and related statistics (`itons()`), for the total number of segregating sites (`segsites()`), and for the tail statistic (`tailstat()`). All these quantities are represented using the `disc_phase_type()` class, whose generator function can also be used for user-tailored discrete phase-type distributions.

The mean and the variance of the any phase-type representation can easily be computed using `mean()` and `var()` respectively. Moreover, a readable summary of a distribution can be obtained using the generic function `summary()`.

This package also contains functions for the density (`dphtype()`), quantiles (`qphtype()`), distribution (`pphtype()`) and random draw generator (`rphtype()`) for both discrete and continuous phase-type distributions. Additionally, the user can calculate and plot the site-frequency spectrum using `sfs()`.

The vignette for the package can be consulted [here](https://aumath-advancedr2019.github.io/ticphasetype/articles/ticphasetype.html). For a quick overview of all included functions, please visit [this link](https://aumath-advancedr2019.github.io/ticphasetype/reference/index.html) . 

## Installation

If you install devtools in your R environment `install.packages("devtools")`, our package will install with the following command.

```
devtools::install_github("aumath-advancedr2019/ticphasetype")
```
Please be aware that our package depends on the [R partitions package](https://cran.r-project.org/web/packages/partitions/index.html), which is only supported by R version >= 3.6.0.
