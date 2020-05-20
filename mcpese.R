#' Accept-reject and Metropolis-Hastings algorithm for Graphical lasso (Glasso) tuning parameter posterior
#'
#' Simulates the posterior distribution of the Glasso tuning parameter. It can be applied for regularization selection of sparse undirected network estimates. This version is compatible with huge package (last tested with v 1.3.2).
#' @param est An object with S3 class "huge".
#' @param n The sample size.
#' @param M The number of accept-reject samples with method = "A-R". The default is 5000.
#' @param g The candidate density in a function form g = function(x)...
#' @param method The sampling method: Accept-Reject "A-R" or Metropolis-Hastings "M-H". The default is "A-R".
#' @param prior The prior distribution with two options: "unif" (if method = "A-R" or "M-H") or "gamma" (only if method = "M-H").
#' @param MH.sampling The style of the Metropolis-Hasting algorithm with two options: "unif" or "random.walk" (default).
#' @param nBurning The number of burn-in iterations. Only applicable when method = "M-H". The default is 1000.
#' @param nSteps The number of iterations. Only applicable when method = "M-H". The default is 5000.
#' @param rhoPriora Shrinkage hyperparameter of the tuning parameter (gamma distribution shape). Only applicable when prior = "gamma". The default is 1.
#' @param rhoPriorb Shrinkage hyperparameter of the tuning parameter (gamma distribution scale). Only applicable when prior = "gamma". The default is 1/10.
#' @param delta the step length of the M-H algorithm. The default is 2.
#' @return A list containing the following components:
#' \itemize{
#' \item est$indx - The vector of indices of the selected tuningparameter values.
#' \item est$rhos - Accepted tuning parameter values.
#' \item est$accept.rate - accept rate.
#' \item est$opt.rho - The smallest of the tuning parameter values which is greater than the expected value of the tuning paramter values sampled either with A-R or M-H.
#' \item est$opt.index - The ordinal of the opt.rho.
#' \item est$n - The sample size.
#' }
#' 
#' @keywords graphical lasso glasso huge network model selection sparse
#' @export
#' @examples
#' library(huge)
#' library(igraph)
#' 
#' set.seed(46023979)
#' 
#' L = huge.generator(d=100, n=120, graph = "hub", g=5)
#' 
#' Y = L$data
#' 
#' nlambda = 50
#' 
#' HugeSolutionPath = huge(Y, method="glasso", nlambda=nlambda)
#' 
#' MCPeSeSelect = mcpese(HugeSolutionPath, n=120)
#' 
#' rhos = MCPeSeSelect$rhos
#' 
#' The function returns the smallest tuning parameter value from the tuning parameter values larger than the simulated mean.
#'
#' optMCPeSelambdaIndx = MCPeSeSelect$opt.index
#'
#' huge.plot(L$theta)
#'
#' title("Ground truth")
#'
#' huge.plot(HugeSolutionPath$path[[optMCPeSelambdaIndx]])
#'
#' title("MCPeSe, accept-rejection")
#'
#'
#' @export
#'
#' @author Markku Kuismin, Mikko J. Sillanpaa
#'
#' @references Kuismin and Sillanpaa (2020) MCPeSe: Monte Carlo penalty selection forgraphical lasso

mcpese = function(est, n=NULL, M=5000, g=NULL, method="A-R", prior="unif", MH.sampling="unif",
                   nBurning=1000, nSteps=5000, rhoPriora=1, rhoPriorb=1/10, delta=2){
  
  # est = huge object estimated with Glasso
  
  # g = candidate density
  
  # prior = prior density
  
  if(is.null(est$loglik)) stop("No log-likelihood")
  
  if(is.null(n)) stop("Set the sample size n")
  
  if(!is.null(g)){
    g = match.fun(g)
    gind = T
  }
  
  if(is.null(g)){
    g = function(rho) 1/(max(rho) - min(rho))
    gind = F
  }
  
  p = ncol(est$data)
  
  rho = est$lambda
  
  nrho = length(rho)
  
  Target = rep(0, nrho)
  
  if(prior == "unif" | method == "A-R"){
    
    for(i in 1:nrho){
      
      Target[i] = p*(p + 1)*log(n*rho[i]/2) - n*rho[i]*sum(abs(est$icov[[i]])) 
      # logarithm of the conditional posterior of the tuning parameter with uniform prior
      
    }
    
    pg = dunif(runif(1, min(rho), max(rho)), min(rho), max(rho), log=T)
    
    Target = Target + 2*pg
    
  }
  
  if(prior == "gamma" & method == "M-H"){
    
    rhoPosta = (rhoPriora + (p*(p + 1)/2))
    
    rhoPostb = rep(0, nrho)
    
    for(i in 1:nrho){
      
      rhoPostb[i] = (rhoPriorb + sum(abs(est$icov[[i]]))/2) 
      
      rhoPostb[i] = n*rhoPostb[i]
      
      Target[i] = stats::dgamma(rho[i], shape=rhoPosta, scale=1/rhoPostb[i], log=T)
      
    }
     
  }
  
  # Accept-reject algorithm:
  
  if(method == "A-R"){
    
    Max = max(Target)*g(rho)
    #Max = max(Target) # (*)
    
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
    #U = runif(M) # (*)
    
    # Hox! U is not the true logartihm transform we need log(exp(Max)*u) but it has the same maximum value and the most extreme
    # upper quantiles are probably close to each other...
    
    Target = Target[Propose]
    
    indx = Propose[U <= Target]
    
    rhos = rho[indx] 
    
    accept.rate = length(indx)/M
    
  }
  
  if(method == "M-H"){
    
    # Initialize rho:
    
    old.state = sample(1:nrho, 1)
    
    rhos = rho[old.state]
    
    indx = old.state
    
    # Metropolis-Hasting sampling:
    
    MCMCSteps = nBurning + nSteps
    
    accept.rate = 0
    
    for(i in 1:MCMCSteps){
      
      # propose new rho:
      
      if(MH.sampling == "unif") new.state = sample(1:nrho, 1) # Jump all over the solution path
      
      if(MH.sampling == "random.walk"){ # Move back and forth of the solution path
        
        new.state = old.state + sample(c(-delta, delta), 1)
        
        if(new.state <= 0) new.state = 1
        
        if(new.state > nrho) new.state = nrho
        
      }
      
      newrho = rho[new.state]
      
      if(prior == "gamma"){
        
        prob = c(0, 
                 stats::dgamma(newrho, shape=rhoPosta, scale=1/rhoPostb[new.state], log=T) -
                   stats::dgamma(rho[old.state], shape=rhoPosta, scale=1/rhoPostb[old.state], log=T)
        )
        
      }
      
      if(prior == "unif"){
        
        prob = c(0, Target[new.state] - Target[old.state])
        
      }
      
      prob = min(prob)
      
      Ind = prob > log(runif(1))
      
      if(Ind){
        rhos = c(rhos, newrho)
        indx = c(indx, new.state)
      }else{
        rhos = c(rhos, rhos[length(rhos)])
        indx = c(indx, old.state)
      }
      
      old.state = new.state
      
      accept.rate = accept.rate + Ind
      
    }
    
    accept.rate = accept.rate/MCMCSteps
    
    rhos = tail(rhos, nSteps)
    
    indx = tail(indx, nSteps)
    
  }
  
  mean.rho = mean(rhos)
  
  opt.rho = min(rho[rho >= mean.rho])
  
  opt.index = which(rho == opt.rho)
  
  est = list()
  
  est$indx = indx
  
  est$rhos = rhos
  
  est$accept.rate = accept.rate
  
  est$opt.rho = opt.rho
  
  est$opt.index = opt.index
  
  est$n = n
  
  return(est)
  
}
