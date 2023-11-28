# Define a function to perform computations for each  window
Pipeline <- function(inSample, outofSample, Update, copulas,
                     Alpha = 0.05, TargetReturn = 0, NumAssets = 8){
  
  # Set seed
  set.seed(2023)
  
  # Setting my pipeline
  Pipe <- map(Update, .f = function(x){
    
    # Create returns matrix
    returns <- inSample[[paste(x)]]
    
    # Fit the GARCH model to the returns data
    fit_garch <- FitGarch(returns)
    
    # Subset the matrix to keep only columns with complete cases
    garch_coef <- Filter(Negate(is.null), fit_garch$garch_coef) # Filtering NULL values 
    unif_dist <- fit_garch$unif_dist
    unif_dist <- unif_dist[, complete.cases(t(unif_dist))] # drop Na columns
    sigma <- fit_garch$sigma
    returns <- returns[, complete.cases(t(sigma))] # drop invalid stocks
    sigma <- sigma[,complete.cases(t(sigma))] # drop Na columns
    
    # Generating Mixture-Copula
    copula_mixture <- tryCatch(
      {
        OptMixtureCopulas(unif_dist, K = 10000, combination = copulas)
      },
      error = function(e) {
        # If an error occurs, adjust uniform dist to have finite limits
        unif_dist <- ifelse(unif_dist < 0.01, 0.01, unif_dist) # avoid convergence issues
        unif_dist <- ifelse(unif_dist > 0.99, 0.99, unif_dist) # avoid convergence issues
        
        # Retry   
        OptMixtureCopulas(unif_dist, K = 10000, combination = copulas)
      }
    )
    
    # Compute simulated standardized residuals using the mixture-copula and GARCH 
    zsim <- ComputeZSim(copula_mixture = copula_mixture, 
                        garch_coef = garch_coef)
    
    # Predict future returns using the GARCH model, simulated residuals, and volatility estimates
    ret_pred <- PredictGarch(returns = returns, 
                             sigma = sigma,
                             zsim = zsim,
                             garch_coef = garch_coef)
    ret_pred <- as.data.frame(ret_pred)
    colnames(ret_pred) <- colnames(returns)
    
    # Perform CVaR optimization to determine the optimal portfolio weights
    weights <- rep(0, ncol(returns))
    names(weights) <- colnames(returns)
    weights <- CVaROptimization(returns = ret_pred,
                                Alpha = Alpha, 
                                TargetReturn = TargetReturn,
                                NumAssets = NumAssets)
    
    # Calculate portfolio returns based on the optimal weights 
    ret_matrix_outofsample <- outofSample[[paste(as.Date(x) + 365)]][,colnames(returns)] # select valid stocks
    portfolio_returns <- as.data.frame(ret_matrix_outofsample  %*%  weights)
  }) %>% 
    bind_rows()
  
  return(Pipe)
}



PipelineParallel <- function(inSample, outofSample, Update, copulas,
                             Alpha = 0.05, TargetReturn = 0, NumAssets = 8){
  
  # Set seed
  set.seed(2023)
  
  # Setting up parallel backend
  registerDoParallel(cores = detectCores())
  
  # Create a foreach loop with parallel execution
  ParallelComputing <- foreach(x = Update, .combine = "bind_rows") %dopar% {
    source("garch_estimate.R")
    source("copula_estimate.R")
    source("portfolio_optimization.R")
    
    # Create returns matrix
    returns <- inSample[[paste(x)]]
    
    # Fit the GARCH model to the returns data
    fit_garch <- FitGarch(returns)
    
    # Subset the matrix to keep only columns with complete cases
    garch_coef <- Filter(Negate(is.null), fit_garch$garch_coef)
    unif_dist <- fit_garch$unif_dist
    unif_dist <- unif_dist[, complete.cases(t(unif_dist))] 
    sigma <- fit_garch$sigma
    returns <- returns[, complete.cases(t(sigma))] 
    sigma <- sigma[,complete.cases(t(sigma))]
    
    # Generating Mixture-Copula
    copula_mixture <- tryCatch(
      {
        OptMixtureCopulas(unif_dist, K = 10000, combination = copulas)
      },
      error = function(e) {
        unif_dist <- ifelse(unif_dist < 0.01, 0.01, unif_dist)
        unif_dist <- ifelse(unif_dist > 0.99, 0.99, unif_dist)
        OptMixtureCopulas(unif_dist, K = 10000, combination = copulas)
      }
    )
    
    # Compute simulated standardized residuals using the mixture-copula and GARCH 
    zsim <- ComputeZSim(copula_mixture = copula_mixture, 
                        garch_coef = garch_coef)
    
    # Predict future returns using the GARCH model, simulated residuals, and volatility estimates
    ret_pred <- PredictGarch(returns = returns, 
                             sigma = sigma,
                             zsim = zsim,
                             garch_coef = garch_coef)
    ret_pred <- as.data.frame(ret_pred)
    colnames(ret_pred) <- colnames(returns)
    
    # Perform CVaR optimization to determine the optimal portfolio weights
    weights <- rep(0, ncol(returns))
    names(weights) <- colnames(returns)
    weights <- CVaROptimization(returns = ret_pred,
                                Alpha = Alpha, 
                                TargetReturn = TargetReturn,
                                NumAssets = NumAssets)
    
    # Calculate portfolio returns based on the optimal weights 
    ret_matrix_outofsample <- outofSample[[paste(as.Date(x) + 365)]][,colnames(returns)] 
    portfolio_returns <- as.data.frame(ret_matrix_outofsample  %*%  weights)
    
    # Retorna o resultado da iteração
    return(data.frame(portfolio_returns))
  }
  
  # Stop parallel backend
  stopImplicitCluster()
  
  return(ParallelComputing)
}