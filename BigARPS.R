
BigARPS = function(Y, M=NULL, nrho=NULL, rho.min.ratio = NULL, g=NULL, verbose=T){
  
  # Accept-Rejest Regularization Selection
  
  # Y = n x p data matrix
  
  # M = nmb of A-R samples (default 1000)
  
  # nrho = nmb of tuning parameters (default 10)
  
  # rho.min.ratio = The smallest value for rho, as a fraction of the uppperbound (MAX) of the 
  # regularization/thresholding parameter which makes all estimates equal to 0. The default value is 0.1
  # See the vigenete of the huge R package
  
  # g = candidate density (default unif(min(rho), max(rho)))
  
  # verbose = Tracking information. Default TRUE
  
  library(glasso)
  
  n = nrow(Y)
  
  p = ncol(Y)
  
  Y = scale(Y)
  
  S = cor(Y)
  
  if(is.null(rho.min.ratio)) rho.min.ratio = 0.1
  
  if(is.null(M)) M = 1000
  
  if(!is.null(g)){
    g = match.fun(g)
    gind = T
  }
  
  if(is.null(g)){
    g = function(rho) 1/(max(rho) - min(rho))
    gind = F
  }
  
  if(is.null(nrho)) nrho = 10
  
  rhomax = max(max(S - diag(p)), -min(S - diag(p)))
  
  rhomin = rho.min.ratio*rhomax
  
  rho = exp(seq(log(rhomax), log(rhomin), length = nrho))
  
  # First compute the solution path with glasso using warm start
  
  Target = rep(0, nrho)
  
  est = glasso(S, rho=rho[1])
  
  Target[1] = p*(p + 1)*log(n*rho[1]/2) - n*rho[1]*sum(abs(est$wi)) 
  
  dir.create("BigARPSSolPath", showWarnings = FALSE)
  
  save(est, file="BigARPSSolPath/est.RDS")
  
  for(i in 2:nrho){
	
    est = glasso(S, rho=rho[i], w.init=est$w, wi.init=est$wi, start="warm")
    
    Target[i] = p*(p + 1)*log(n*rho[i]/2) - n*rho[i]*sum(abs(est$wi)) 
    # logarithm of the conditional posterior of the tuning parameter
    
    saveRDS(est, file=paste("BigARPSSolPath/est", i, ".RDS", sep=""))
    
    if (verbose) {
      
      cat("Running glasso", round(100*i/nrho, 2), "%", "\r")
      
    }
    
  }
  
  # Accept-reject algorithm:
  
  pg = 1/(max(rho) - min(rho))
  
  Target = Target + 2*log(pg)
  
  Max = max(Target)*g(rho)
  
  Propose = sample(1:nrho, M, replace=T)
  
  if(gind == T){
    
    Mg = max(g(rho))
    
    s = rho[Propose]
    
    s = Propose[runif(length(s)) <= g(s)/Mg]
    
    IndProp = sample(1:nrho, 1, replace = T)
    
    while(length(s) < M){
      
      if(runif(1) <= g(rho[IndProp])/Mg) s = c(s, IndProp)
      
      IndProp = sample(1:nrho, 1, replace = T)
      
    }
    
    Propose = s
    
  }
  
  U = runif(M, 0, Max) 
  
  # Hox! U is not the true logartihm transform we need log(exp(Max)*u) but it has the same maximum value and the most extreme
  # upper quantiles are probably close to each other...
  
  Target = Target[Propose]
  
  indx = Propose[U <= Target] # No log-likelihood at all!
  
  Sampledrhos = rho[indx] 
  
  Results = list()
  
  Results$indx = indx
  
  Results$Sampledrhos = Sampledrhos
  
  Results$rhos = rho
  
  Results$accept.rate = length(indx)/M
  
  Results$n = n
  
  return(Results)
  
}
