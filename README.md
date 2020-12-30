SpReg: A Toolkit of Spatial Econometrics with Misaligned Data
================
*[Guillume Pouliot](https://sites.google.com/site/guillaumeallairepouliot) and [Zhen Xie](https://github.com/zhenxie23/)*

Introduction
------------

This R package implements regression analysis for spatially misaligned data, where the geographic locations of the outcome variable and explanatory variables do not coincide, as developed in [Pouliot (2020)](https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnxndWlsbGF1bWVhbGxhaXJlcG91bGlvdHxneDoxN2QzNjYwNmQ5ODczYjE). We implement estimation and inference for three complementary estimators: a two-step bootstrap estimator, a minimum-distance estimator, and a one-step Quasi Maximum Likelihood estimator. These three methods trade off statistical power for weaker assumptions on the covariance structure of regression residuals.

Installation
-----------------------------

`SpReg` can be installed from its GitHub repository via
```r
devtools::install_github("zhenxie23/SpReg")
```

Replication Example: [Madsen et al. (2008)](https://onlinelibrary.wiley.com/doi/abs/10.1002/env.888)
------------------------------
We use the dataet of [Madsen et al. (2008)](https://onlinelibrary.wiley.com/doi/abs/10.1002/env.888). There are 558 observations over an area of 400,000 squared kilometers.
```r
library(SpReg)
rivers <- read.csv(system.file("extdata", "rivers.csv", package = "SpReg"))
rivers <- rivers[which(rivers$FOR_NLCD<100),]
```
The simulation is implemented as follows: for each run, the 558 locations in the data set are randomly spilt into two equal halves. The covariate R = logit(for.nlcd) is assumed to have been observed on one half, and the out come variable Y = log(cl) is assumed to have been observe on the other. This creates a misaligned data set. For each round of the cross-validation exercise, all estimates, along with their standar errors, are calculated. The performance of estimators can then be compare with each other.
```r
train_set <- sample(1:558,277);
test_set <- setdiff(1:558,train_set);
rivers_train <- rivers[train_set, ];
rivers_test <- rivers[test_set, ];

DatR <- rivers_train[,c("FOR_NLCD", "LAT_DD", "LON_DD")];
DatR$X <- log((DatR$FOR_NLCD)/(100-DatR$FOR_NLCD));
sp::coordinates(DatR) <- ~LON_DD+LAT_DD;
sp::proj4string(DatR) =  "+proj=longlat +datum=WGS84";

DatY <- rivers_test[, c("CL","LAT_DD","LON_DD")];
DatY$Y <- log(DatY$CL);
sp::coordinates(DatY) <- ~LON_DD+LAT_DD;
sp::proj4string(DatY) =  "+proj=longlat +datum=WGS84";
```

*Estimation and Inference for Two-Step Estimators*
```r
TwoStep_Results <- TwoStepBootstrap(DatR, "X", DatY, "Y", "Exp", FALSE, FALSE, cutoff = 295, cutoff_u = 40);
```
