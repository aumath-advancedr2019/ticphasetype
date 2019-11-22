# ticphasetype [![Build Status](https://travis-ci.com/aumath-advancedr2019/ticphasetype.svg?branch=master)](https://travis-ci.com/aumath-advancedr2019/ticphasetype)

## Outline for package

### Density, cummulative probability, quantile function and random generator for continuous and discrete phase-type distributions:
* dphasetype()
* pphasetype()
* qphasetype()
* rphasetype()

### Phase-Type matrices for basic Kingsman Coalescent 
* T_mrca 
* T_total 
* Kingsman

### Phase-type representation of different summary statistics 
* mean and variance
* S_total 
* i-tons 
* tail statistic 
* sfs (site frequency spectrum)

## Roles

* pm: Iker.
* de: Carl.
* qe: Tomas. 





## Installation

If you install devtools in your R environment `install.packages("devtools")`, our package will install with the following command.

```
devtools::install_github("aumath-advancedr2019/ticphasetype")
```
Please be aware that our package depends on the [R partitions package](https://cran.r-project.org/web/packages/partitions/index.html), which is only supported by R version >= 3.6.0.
