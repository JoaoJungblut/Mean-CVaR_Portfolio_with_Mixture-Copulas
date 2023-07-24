
# Rolling window estimation function
RollingWindowEstimation <- function(returns, We = 252, Wt = ncol(returns), K = 10000) {
  
  # Initialize empty data frames to store the results
  names_vector <- names(returns)[-1]   # Asset names for reference
  weights <- matrix(nrow = Wt, ncol = N) # Create a matrix to store the weights for each asset in the portfolio
  colnames(weights) <- names_vector # Set the column names of the weights matrix as the asset names
  weights[1:Wt,] <-  0 # Initialize the first K rows of the weights matrix as all zeros
  portfolio_returns <- matrix(nrow = Wt, ncol = 1)  # Matrix to store portfolio returns
  portfolio_returns[1:Wt, ] <- 0  # Initialize the first K rows as zero
  weights_benchmark <- matrix(nrow = Wt, ncol = N) # Create a matrix to store the weights for each asset in the Gaussian benchmark portfolio
  colnames(weights_benchmark) <- names_vector # Set the column names of the weights matrix as the asset names
  weights_benchmark[1:Wt,] <-  0 # Initialize the first K rows of the weights matrix as all zeros
  benchmark_portfolio_returns <- matrix(nrow = Wt, ncol = 1)  # Matrix to store portfolio returns
  benchmark_portfolio_returns[1:Wt, ] <- 0  # Initialize the first K rows as zero
  
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
    
    # Optimize the mixture of copulas using the uniform distribution from the GARCH model
    copulas_mixture <- OptMixtureCopulas(unif_dist = fit_garch$unif_dist, K = K)
    
    # Compute simulated standardized residuals using the optimized copula mixture and GARCH coefficients
    zsim <- ComputeZSim(copula_mixture = copulas_mixture, garch_coef = fit_garch$garch_coef)
    
    # Generate Gaussian Copula for benchmark portfolio
    copulas_gauss <- GaussCopula(unif_dist = fit_garch$unif_dist, K = K)
    
    # Compute simulated standardized residuals from Gaussian Copula
    gsim <- ComputeZSim(copula_mixture = copulas_gauss, garch_coef = fit_garch$garch_coef)
    
    # Predict future returns using the GARCH model, simulated residuals, and volatility estimates
    ret_pred <- PredictGarch(returns = ret_matrix_insample, 
                             sigma = fit_garch$sigma,
                             zsim = zsim,
                             garch_coef = fit_garch$garch_coef)
    
    ret_benchmark_pred <- PredictGarch(returns = ret_matrix_insample, 
                                       sigma = fit_garch$sigma,
                                       zsim = gsim,
                                       garch_coef = fit_garch$garch_coef)
    
    # Perform CVaR optimization to determine the optimal portfolio weights
    weights[i, names_vector[assets_with_valid_returns]] <- CVaROptimization(returns = ret_pred)
    weights[i, names_vector[!assets_with_valid_returns]] <- 0
    
    weights_benchmark[i, names_vector[assets_with_valid_returns]] <- CVaROptimization(returns = ret_benchmark_pred)
    weights_benchmark[i, names_vector[!assets_with_valid_returns]] <- 0
    
    # Convert the realized returns data to a matrix format
    ret_matrix_outofsample <- as.matrix(returns[i, -1])
    ret_matrix_outofsample[, names_vector[!assets_with_valid_returns]] <- 0
    
    # Calculate the portfolio returns based on the optimal weights
    portfolio_returns[i, ] <- RetPortfolio(returns = ret_matrix_outofsample,  
                                           weights = rbind(weights[i, ])) - 0.0003 # minus the transaction costs
    
    benchmark_portfolio_returns[i, ] <- RetPortfolio(returns = ret_matrix_outofsample,  
                                                     weights = rbind(weights_benchmark[i, ])) - 0.0003 # minus the transaction cost
  }
  
  # Apply the rolling window computation to each window
  purrr::map((We + 1):Wt, process_window)
  
  return(list(
    portfolio_returns = portfolio_returns,
    benchmark_portfolio_returns = benchmark_portfolio_returns,
  ))
}

# Sample usage:
result <- RollingWindowEstimation(returns, We=252, Wt=263, K=1000, Mixture = TRUE)


# Rolling window estimation function
RollingWindowEstimation <- function(returns, 
                                    We = 252,
                                    Wt = ncol(returns), 
                                    K = 10000,
                                    Mixture = TRUE) {
  
  if(Mixture == TRUE){
    # Initialize empty data frames to store the results
    names_vector <- names(returns)[-1]   # Asset names for reference
    weights <- matrix(nrow = Wt, ncol = N) # Create a matrix to store the weights for each asset in the portfolio
    colnames(weights) <- names_vector # Set the column names of the weights matrix as the asset names
    weights[1:We,] <-  0 # Initialize the first K rows of the weights matrix as all zeros
    portfolio_returns <- matrix(nrow = Wt, ncol = 1)  # Matrix to store portfolio returns
    portfolio_returns[1:We, ] <- 0  # Initialize the first K rows as zero
    
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
      
      # Optimize the mixture of copulas using the uniform distribution from the GARCH model
      copulas_mixture <- OptMixtureCopulas(unif_dist = fit_garch$unif_dist, K = K)
      
      # Compute simulated standardized residuals using the optimized copula mixture and GARCH coefficients
      zsim <- ComputeZSim(copula_mixture = copulas_mixture, garch_coef = fit_garch$garch_coef)
      
      # Predict future returns using the GARCH model, simulated residuals, and volatility estimates
      ret_pred <- PredictGarch(returns = ret_matrix_insample, 
                               sigma = fit_garch$sigma,
                               zsim = zsim,
                               garch_coef = fit_garch$garch_coef)
      
      # Perform CVaR optimization to determine the optimal portfolio weights
      weights[i, names_vector[assets_with_valid_returns]] <- CVaROptimization(returns = ret_pred)
      weights[i, names_vector[!assets_with_valid_returns]] <- 0
    
      # Convert the realized returns data to a matrix format
      ret_matrix_outofsample <- as.matrix(returns[i, -1])
      ret_matrix_outofsample[, names_vector[!assets_with_valid_returns]] <- 0
      
      # Calculate the portfolio returns based on the optimal weights
      portfolio_returns[i, ] <- RetPortfolio(returns = ret_matrix_outofsample,  
                                             weights = rbind(weights[i, ])) - 0.0003 # minus the transaction costs
      
    }
    
    # Apply the rolling window computation to each window
    purrr::map((We + 1):Wt, process_window)
    result <- data.frame(date = returns[1:Wt,]$date, 
                         portfolio_return = portfolio_returns[1:Wt,])
    
  } else{
    # Initialize empty data frames to store the results
    names_vector <- names(returns)[-1]   # Asset names for reference
    weights_benchmark <- matrix(nrow = Wt, ncol = N) # Create a matrix to store the weights for each asset in the Gaussian benchmark portfolio
    colnames(weights_benchmark) <- names_vector # Set the column names of the weights matrix as the asset names
    weights_benchmark[1:We,] <-  0 # Initialize the first K rows of the weights matrix as all zeros
    benchmark_portfolio_returns <- matrix(nrow = Wt, ncol = 1)  # Matrix to store portfolio returns
    benchmark_portfolio_returns[1:We, ] <- 0  # Initialize the first K rows as zero
    
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
      
      # Generate Gaussian Copula for benchmark portfolio
      copulas_gauss <- GaussCopula(unif_dist = fit_garch$unif_dist, K = K)
      
      # Compute simulated standardized residuals from Gaussian Copula
      gsim <- ComputeZSim(copula_mixture = copulas_gauss, garch_coef = fit_garch$garch_coef)
      
      # Predict future returns using the GARCH model, simulated residuals, and volatility estimates
      ret_benchmark_pred <- PredictGarch(returns = ret_matrix_insample, 
                                         sigma = fit_garch$sigma,
                                         zsim = gsim,
                                         garch_coef = fit_garch$garch_coef)
      
      # Perform CVaR optimization to determine the optimal portfolio weights
      weights_benchmark[i, names_vector[assets_with_valid_returns]] <- CVaROptimization(returns = ret_benchmark_pred)
      weights_benchmark[i, names_vector[!assets_with_valid_returns]] <- 0
      
      # Convert the realized returns data to a matrix format
      ret_matrix_outofsample <- as.matrix(returns[i, -1])
      ret_matrix_outofsample[, names_vector[!assets_with_valid_returns]] <- 0
      
      # Calculate the portfolio returns based on the optimal weights
      benchmark_portfolio_returns[i, ] <- RetPortfolio(returns = ret_matrix_outofsample,  
                                                       weights = rbind(weights_benchmark[i, ])) - 0.0003 # minus the transaction cost
    }
    
    # Apply the rolling window computation to each window
    purrr::map((We + 1):Wt, process_window)
    result <- data.frame(date = returns[1:Wt,]$date, 
                         portfolio_return = benchmark_portfolio_returns[1:Wt,])
    
  }
  
  return(result)
}
