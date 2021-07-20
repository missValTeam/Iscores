# Iscores: scoring imputations methods


## Overview

Iscores is a package intended to provide a framework for scoring imputations methods. It implements the scores described in [Michel, Naef, Spohn and Meinshausen. 2021](https://arxiv.org/abs/2106.03742) . Examples of use of the library are shown below.


## Installation

The package should be (soon) available on CRAN, To install the package from github you can run

``` r
install.packages("devtools")
devtools::install_github("missValTeam/Iscores")
```

## Examples: 


```r
n <- 20
X <- cbind(rnorm(n),rnorm(n))
X.NA <- X
X.NA[,1] <- ifelse(runif(n)<=0.2, NA, X[,1])

imputations <- list()

imputations[[1]] <- lapply(1:5, function(i) {
  X.loc <- X.NA
  X.loc[is.na(X.NA[,1]),1] <- mean(X.NA[,1],na.rm=TRUE)
  return(X.loc)
})

imputations[[2]] <- lapply(1:5, function(i) {
  X.loc <- X.NA
  X.loc[is.na(X.NA[,1]),1] <- sample(X.NA[!is.na(X.NA[,1]),1], size = sum(is.na(X.NA[,1])), replace = TRUE)
  return(X.loc)
})

methods <- c("mean","sample")

Iscores(imputations = imputations, methods = methods, X.NA = X.NA)
```


## Issues

To report an issue, please use the [issue tracker](https://github.com/missValTeam/Iscores/issues) on github.com.
