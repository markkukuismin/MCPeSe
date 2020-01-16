
ARPS = function(est, n=NULL, M=NULL, g=NULL){
  
  # Accept-Rejest Regularization Selection
  
  # est = huge object estimated with Glasso
  
  # n = sample size
  
  # M = nmb of A-R samples
  
  # g = candidate density
  
  if(is.null(est$loglik)) stop("No log-likelihood")
  
  if(is.null(n)) stop("Set the sample size n")
  
  if(is.null(M)) M = 1000
  
  if(!is.null(g)){
    g = match.fun(g)
    gind = T
  }
  
  if(is.null(g)){
    g = function(rho) 1/(max(rho) - min(rho))
    gind = F
  }
  
  p = ncol(est$data)
  
  RemoveInfinite = c(which(est$loglik == Inf), which(est$loglik == -Inf)) 
  
  rho = est$lambda
  
  nlambda = length(rho)
  
  if(length(RemoveInfinite) > 0){
    
    rho = rho[-RemoveInfinite]
    
    nlambda = length(rho)
    est$path = est$path[-RemoveInfinite]
    est$lambda = est$lambda[-RemoveInfinite]
    est$icov = est$icov[-RemoveInfinite]
    est$df = est$df[-RemoveInfinite]
    est$sparsity = est$sparsity[-RemoveInfinite]
    est$loglik = est$loglik[-RemoveInfinite]
    
  }
  
  # Accept-reject algorithm:
  
  Target = rep(0,nlambda)
  
  for(i in 1:nlambda){
    
    Target[i] = p*(p + 1)*log(n*est$lambda[i]/2) - n*est$lambda[i]*sum(abs(est$icov[[i]])) 
    # logarithm of the conditional posterior of the tuning parameter
    
  }
  
  pg = 1/(max(rho) - min(rho))
  
  Target = Target + 2*log(pg)
  
  Max = max(Target)*g(rho)
  
  Propose = sample(1:nlambda, M, replace=T)
  
  if(gind == T){
  
	Mg = max(g(rho))

	s = rho[Propose]

	s = Propose[runif(length(s)) <= g(s)/Mg]

	IndProp = sample(1:nlambda, 1, replace = T)

	while(length(s) < M){
  
		if(runif(1) <= g(rho[IndProp])/Mg) s = c(s, IndProp)
  
		IndProp = sample(1:nlambda, 1, replace = T)
  
	}
  
	Propose = s
  
  }
  
  U = runif(M, 0, Max) 
  
  # Hox! U is not the true logartihm transform we need log(exp(Max)*u) but it has the same maximum value and the most extreme
  # upper quantiles are probably close to each other...
  
  Target = Target[Propose]
  
  indx = Propose[U <= Target] # No log-likelihood at all!
  
  rhos = rho[indx] 
  
  est = list()
  
  est$indx = indx
  
  est$rhos = rhos
  
  est$accept.rate = length(indx)/M
  
  est$n = n
  
  return(est)
  
}
