ComputePerformance <- function(returns) {
  
  # ComputePerformance: Function to measure portfolio performance using various performance ratios.
  # Inputs:
  #   returns: A data frame containing asset returns with a 'date' column and individual asset columns.
  # Output:
  #   A list of computed performance ratios.
  
  
  # Compute annualized mean excess return
  ret_annualized <- PerformanceAnalytics::Return.annualized(returns)
  
  # Compute annualized standard deviation
  std_dev_annualized <- PerformanceAnalytics::StdDev.annualized(returns)
  
  # Compute annualized Sharpe Ratio
  Sharpe_annualized <- PerformanceAnalytics::SharpeRatio.annualized(returns)
  
  # Compute Sortino Ratio
  Sortino <- PerformanceAnalytics::SortinoRatio(returns)
  
  # Compute Omega Sharpe Ratio
  omega <- PerformanceAnalytics::OmegaSharpeRatio(returns)
  
  # Compute VaR (Value at Risk) at 97.5% confidence level
  var_975 <- PerformanceAnalytics::VaR(returns, p = 0.975, method = "historical")
  
  # Compute Conditional VaR (CVaR) at 97.5% confidence level
  cvar_975 <- PerformanceAnalytics::ES(returns, p = 0.975, method = "historical")
  
  # Compute semi-deviation
  semi_deviation <- PerformanceAnalytics::SemiDeviation(returns)
  
  # Compute worst drawdown
  worst_drawdown <- PerformanceAnalytics::maxDrawdown(returns)
  
  # Create a list of performance ratios
  performance_ratios <- list(ret_annualized,
                             std_dev_annualized,
                             Sharpe_annualized,
                             Sortino,
                             omega,
                             var_975,
                             cvar_975,
                             semi_deviation,
                             worst_drawdown)
  
  return(performance_ratios)
}
