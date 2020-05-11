![GitHub](https://img.shields.io/github/license/markkukuismin/ARPS)

# ARPS
Accept-Reject Penlaty Selection (ARPS) for graphical lasso.

# Example

```r
library(huge)

source("ARPS_v2.R")

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

HugeData = huge.generator(n=n, d=p, graph=Model) # Just the precision matrix corresponding to the graphical model of
# interest is needed. Data is simulated later.

nlambda = 100

##########################################################

# Compute Glasso solution path:

L = huge(HugeData$data, nlambda=nlambda, method="glasso")

##########################################################

# Run the A-R selection (uniform prior),

ARSelect = ARPS_v2(L, n=n, M=1000)

names(ARSelect)

rhos = ARSelect$rhos

ARSelect$accept.rate

##########################################################

# Run the M-H selection (uniform prior),

MHSelect = ARPS_v2(L, n=n, nSteps=1000, method="M-H")

names(MHSelect)

rhos = MHSelect$rhos

MHSelect$accept.rate

plot(rhos, type="l")

##########################################################
```

![MHMCMCSample](https://user-images.githubusercontent.com/40263834/81556858-b38c1880-9393-11ea-8eda-50274a57d9a6.png)

```r
# Either use the mean of rho values...

mean(rhos)

#ThetaARSelect = huge(Y, lambda=mean(rhos), method="glasso")
#ThetaAR = as.matrix(ThetaARSelect$icov[[1]])

# ... or pick the smallest tuning parameter value from the solution path
# which is larger or equal to the mean value:

optARrhoIndx = ARSelect$opt.index # This is sup{i : rho[i] >= mean(rhos)}

optMHrhoIndx = MHSelect$opt.index

huge.plot(HugeData$theta)

title("Ground truth")

##########################################################
```

![GroundTruthGraph](https://user-images.githubusercontent.com/40263834/81557077-0f56a180-9394-11ea-8827-84b046cbf667.png)

```r
huge.plot(L$path[[optARrhoIndx]])

title("Accept-Reject sampling")

##########################################################
```
![ARGraph](https://user-images.githubusercontent.com/40263834/81557112-1d0c2700-9394-11ea-9883-a7e72507c12f.png)

```r
huge.plot(L$path[[optMHrhoIndx]])

title("Metropolis-Hastings sampling")
```
![MHGraph](https://user-images.githubusercontent.com/40263834/81557140-272e2580-9394-11ea-8cb0-054910e678b2.png)

# Reference

The ARPS method is described in:

Kuismin and Sillanpaa (manuscript) "ARPS: Accept-reject penalty selection for graphical lasso".

File "CodeCollection.zip" is a collection of scripts used to prepare the material in this paper.

