RollingWindowEstimation <- function(returns, 
                                    We = 252,
                                    Wt = ncol(returns), 
                                    K = 10000,
                                    Mixture = TRUE) {
  
  # RollingWindowEstimation: Function to perform rolling window estimation of optimal portfolio weights
  # Inputs:
  #   returns: A data frame containing returns data with the 'date' column and each asset's returns as other columns.
  #   We: Window size (default = 252) representing the number of periods for each rolling window.
  #   Wt: Total number of columns in the 'returns' data frame (default is the number of columns in 'returns').
  #   K: Number of iterations for copula fitting (default = 10000).
  #   Mixture: Boolean flag indicating whether to use a mixture copula (default = TRUE).
  # Output:
  #   A data frame with 'date' and 'portfolio_return' columns, containing the portfolio returns for each rolling window.
  
  
  set.seed(64)
  
  
  # Define a function to perform computations for each rolling window
  process_window <- function(i) {
    t1 <- i - We
    t2 <- i - 1
    
    # Convert the in-sample returns data to a matrix format
    ret_matrix_insample <- as.matrix(returns[t1:t2, -1])
    
    # Create a logical vector indicating if each asset has sufficient data
    assets_with_valid_returns <- !colMeans(is.na(ret_matrix_insample))
    
    # Subset the returns matrix and asset names based on assets with sufficient data
    ret_matrix_insample <- ret_matrix_insample[, assets_with_valid_returns]
    
    # Fit the GARCH model to the returns data
    fit_garch <- FitGarch(returns = ret_matrix_insample)
    
    # Choose the appropriate copula based on Mixture
    association_measure <- if (Mixture) {
      OptMixtureCopulas(unif_dist = fit_garch$unif_dist, K = K)
    } else {
      GaussCopula(unif_dist = fit_garch$unif_dist, K = K)
    }
    
    # Compute simulated standardized residuals using the optimized copula mixture and GARCH coefficients
    zsim <- ComputeZSim(copula_mixture = association_measure, 
                        garch_coef = fit_garch$garch_coef)
    
    # Predict future returns using the GARCH model, simulated residuals, and volatility estimates
    ret_pred <- PredictGarch(returns = ret_matrix_insample, 
                             sigma = fit_garch$sigma,
                             zsim = zsim,
                             garch_coef = fit_garch$garch_coef)
    
    # Perform CVaR optimization to determine the optimal portfolio weights
    weights <- CVaROptimization(returns = ret_pred)
    
    return(weights)
  }
  
  # Apply the rolling window computation to each window and store results in a list
  all_weights <- purrr::map((We + 1):Wt, process_window)
  
  # Combine the results into a single data frame
  weights <- do.call(rbind, all_weights)
  
  # Create a matrix with realized returns for out-of-sample period
  ret_matrix_outofsample <- as.matrix(returns[(We + 1):Wt, -1])
  
  # Calculate portfolio returns based on the optimal weights for each rolling window
  portfolio_returns <- apply(ret_matrix_outofsample, 1, function(row) {
    # Remove rows with NAs from ret_matrix_outofsample
    row_without_na <- row[complete.cases(row)]
    # Calculate portfolio return for this window
    portfolio_return <- sum(row_without_na * weights)
    return(portfolio_return)
  }) - 0.0003 # minus the transaction costs
  
  # Combine dates with portfolio returns into a data frame
  result <- data.frame(date = returns[(We + 1):Wt, "date"], 
                       portfolio_return = portfolio_returns)
  
  return(result)
}

