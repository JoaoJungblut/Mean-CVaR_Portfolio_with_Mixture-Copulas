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
    
    message(paste(i - We, "of", Wt - We))
    
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
    
    # Define a function to compute the association measure with error handling
    association_measure_with_error_handling <- function(unif_dist, Mixture, K) {
      tryCatch(
        {
          if (Mixture) {
            association_measure <- OptMixtureCopulas(unif_dist = unif_dist, K = K)
          } else {
            association_measure <- GaussCopula(unif_dist = unif_dist, K = K)
          }
        },
        error = function(e) {
          # In case of an error, return NULL
          warning("An error occurred while computing the association measure:", 
                  conditionMessage(e))
          association_measure <- NULL 
        }
      )
      
      return(association_measure)
    }
    
    # Compute the association measure with error handling
    association_measure <- association_measure_with_error_handling(unif_dist = fit_garch$unif_dist,
                                                                   Mixture = Mixture, 
                                                                   K = K)
    
    # Check if association_measure is NULL (i.e., an error occurred)
    if (is.null(association_measure)) {
      # Set weights to zero if an error occurred
      weights <- rep(0, ncol(returns[,-1]))
    } else {
      # Compute simulated standardized residuals using the optimized copula mixture and GARCH coefficients
      zsim <- ComputeZSim(copula_mixture = association_measure, 
                          garch_coef = fit_garch$garch_coef)
      
      # Predict future returns using the GARCH model, simulated residuals, and volatility estimates
      ret_pred <- PredictGarch(returns = ret_matrix_insample, 
                               sigma = fit_garch$sigma,
                               zsim = zsim,
                               garch_coef = fit_garch$garch_coef)
      
      # Perform CVaR optimization to determine the optimal portfolio weights
      names_vector <- names(returns[,-1])
      weights <- matrix(nrow = 1, ncol = ncol(returns[,-1])) 
      colnames(weights) <- names_vector
      weights[1, names_vector[assets_with_valid_returns]] <- CVaROptimization(returns = ret_pred)
      weights[1, names_vector[!assets_with_valid_returns]] <- 0
    }
    
    return(weights)
  }
  
  # Apply the rolling window computation to each window and store results in a list
  all_weights <- purrr::map((We + 1):Wt, process_window)
  
  # Combine the results into a single data frame
  weights <- do.call(rbind, all_weights)
  
  # Create a matrix with realized returns for out-of-sample period
  ret_matrix_outofsample <- as.matrix(returns[(We + 1):Wt, -1])
  
  # Replace NA values with 0
  ret_matrix_outofsample[is.na(ret_matrix_outofsample)] <- 0
  
  # Calculate portfolio returns based on the optimal weights 
  portfolio_returns <- rowSums(ret_matrix_outofsample * weights) - 0.0003 # minus the transaction costs
  
  # Combine dates with portfolio returns into a data frame
  result <- data.frame(date = returns[(We + 1):Wt, "date"], 
                       portfolio_return = portfolio_returns)
  
  return(result)
}

