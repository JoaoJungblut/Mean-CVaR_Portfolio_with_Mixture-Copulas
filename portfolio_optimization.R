CVaROptimization <- function(returns, Alpha = 0.025, TargetReturn = 0.0003) { # Target return must consider transaction costs
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


NaiveDiversification <- function(returns) {
  # Calculate the number of assets
  num_assets <- ncol(returns) - 1  # Subtract 1 for the 'date' column
  
  # Calculate the equal weights for each asset
  weight_per_asset <- 1 / num_assets
  
  # Extract the returns data without the 'date' column
  returns_df <- returns[, -1]  # Drop the 'date' column
  
  # Calculate portfolio returns for each date
  portfolio_returns <- rowSums(returns_df, na.rm = TRUE) * weight_per_asset
  
  # Create a new data frame with the date and portfolio returns
  portfolio_returns_df <- data.frame(date = returns$date, portfolio_return = portfolio_returns)
  
  return(portfolio_returns_df)
}

