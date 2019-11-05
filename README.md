# ARRSe
Accept-Reject Regularization Selection (ARRSe) for graphical lasso.

# Example

ARRSe is compatible with the R package "huge".

```r

library(huge)

source("ARRSe.R")

##########################################################

# Initialize seed number:

#seed = Sys.time()

#seed = as.integer(seed)

#seed = seed %% 100000

seed = 59040

set.seed(seed)

##########################################################

# Graphical model simulation:

p = 200

n = 210

Model = "hub"

HugeData = huge.generator(n=n, d=p, graph=Model)

nlambda = 100

##########################################################

# Compute Glasso solution path:

L = huge(HugeData$data, nlambda=nlambda, method="glasso")

##########################################################

# Run ARRSe:

ARSelect = ARRSe(L, n=n, M=1000)

names(ARSelect)

rhos = ARSelect$rhos

ARSelect$accept.rate

# Either use the mean of rho values...

mean(rhos)

#ThetaARSelect = huge(Y, lambda=mean(rhos), method="glasso")
#ThetaARRSe = as.matrix(ThetaARSelect$icov[[1]])

# ... or pick the value from the solution path which is closest to the mean value:

d = abs(mean(rhos) - L$lambda)

optARRSelambdaIndx = which.min(d)[length(which.min(d))] # This is suppremum; lambda sequence is decreasing

huge.plot(HugeData$theta)

title("Ground truth")
```
![GitHubDemoGroundTruth](https://user-images.githubusercontent.com/40263834/68211077-87a46a80-ffdf-11e9-915b-fef820af900e.png)

```r
huge.plot(L$path[[optARRSelambdaIndx]])

title("ARRSe")
```
![GitHubDemoARRSe](https://user-images.githubusercontent.com/40263834/68211099-95f28680-ffdf-11e9-95e1-9e5514f871e6.png)

