CVaROptimization <- function(returns, Alpha = 0.025, TargetReturn = 0) {
  frontierSpec <- fPortfolio::portfolioSpec()  # Portfolio specification for optimization
  fPortfolio::setType(frontierSpec) <- "CVaR"  # Set portfolio type as CVaR
  fPortfolio::setSolver(frontierSpec) <- "solveRglpk.CVAR"  # Use linear programming solver for CVaR optimization
  fPortfolio::setAlpha(frontierSpec) <- Alpha  # Set CVaR alpha level as 0.025 (CVaR_0.975)
  fPortfolio::setTargetReturn(frontierSpec) <- TargetReturn  # Set the daily target return constraint
  
  # Create time series object
  returnfPort <- as.timeSeries(returns)
  
  # Optimizing portfolio using K simulated returns for each asset
  frontier1g <- fPortfolio::efficientPortfolio(data = returnfPort,
                                               spec = frontierSpec,
                                               constraints = "LongOnly")
  
  # Storing resulting weights
  cvar_opt <- rbind(fPortfolio::getWeights(frontier1g))
  
  # Return portfolio weights
  return(cvar_opt)
}


RetPortfolio <- function(returns, weights){
  
  portfolio_returns <- returns %*% t(weights)
  
  return(portfolio_returns)
}



