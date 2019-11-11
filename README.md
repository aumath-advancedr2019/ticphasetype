# ticphasetype [![Build Status](https://travis-ci.com/aumath-advancedr2019/ticphasetype.svg?branch=master)](https://travis-ci.com/aumath-advancedr2019/ticphasetype)

## Outline for package

Define the density, cummulative probability and quantile functions for continuous and discrete phase-type distributions:
* dphasetype()
* pphasetype()
* qphasetype()
* rphasetype()

Create phase-type representation of different summary statistics (T_mrca, T_total, S_total, n-tons...) for any number of n. 

Additional analyses:
* Theta estimators (mutation rate)
* (maybe) We will try to implement summary statistics (e.g. the 15-lumped tail v.s. normalized singletons, which together with LDA can differentiate the β-coalescent from K+Exp-coalescent).

## Roles

* pm: Iker.
* de: Carl.
* qe: Tomas. 





## Installation

If you install devtools in your R environmen install.packages("devtools"), you will be able to install our package with the following command.
`devtools::install_github("aumath-advancedr2019/ticphasetype")`
