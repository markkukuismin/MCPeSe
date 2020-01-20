![GitHub](https://img.shields.io/github/license/markkukuismin/ARPS)

# ARPS
Accept-Reject Penlaty Selection (ARPS) for graphical lasso.

# Example

ARPS is compatible with the R package "huge".

```r

library(huge)

source("ARPS.R")

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

# Run ARPS:

# The default proposal distribution is 1/(max(rho) - min(rho)). 

# Other density functions can also be applied: g = function(x) ...

ARSelect = ARPS(L, n=n, M=1000)

names(ARSelect)

[1] "indx"        "rhos"        "accept.rate" "n"  

rhos = ARSelect$rhos

ARSelect$accept.rate

[1] 0.341

# Either use the mean of rho values...

mean(rhos)

[1] 0.2397052

#ThetaARSelect = huge(Y, lambda=mean(rhos), method="glasso")
#ThetaARPS = as.matrix(ThetaARSelect$icov[[1]])

# ... or pick the value from the solution path which is closest to the mean value:

d = abs(mean(rhos) - L$lambda)

optARPSlambdaIndx = which.min(d)[length(which.min(d))]

huge.plot(HugeData$theta)

title("Ground truth")
```
![GitHubDemoGroundTruth](https://user-images.githubusercontent.com/40263834/68211077-87a46a80-ffdf-11e9-915b-fef820af900e.png)

```r
huge.plot(L$path[[optARPSlambdaIndx]])

title("ARPS")
```
![GitHubDemoARPS](https://user-images.githubusercontent.com/40263834/72513292-9c11a880-3855-11ea-9608-269a259fcd41.jpeg)

# Reference

The ARPS method is described in:

Kuismin and Sillanpaa (manuscript) "ARPS: Accept-reject penalty selection for graphical lasso".

File "CodeCollection.zip" is a collection of scripts used to prepare the material in this paper.

