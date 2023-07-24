FitGarch <- function(returns){
  
  # Initialize the GARCH specification
  mod_garch <- try(rugarch::ugarchspec(variance.model = list(model = "sGARCH", 
                                                             garchOrder = c(1, 1),
                                                             variance.targeting = TRUE),
                                       mean.model = list(armaOrder = c(1, 0)),
                                       distribution.model = "sstd"),
                   silent = TRUE)
  
  # Matrix to save model coefficients, standard deviation and residuals
  unif_dist <- garch_pred <- sigma <- residuals <- matrix(nrow = nrow(returns), 
                                                          ncol = ncol(returns))
  garch_coef <- vector("list", length = ncol(returns))
  
  # Loop over each asset
  for(j in 1:ncol(returns)){
    # Fitting GARCH model
    garch_fit <- try(ugarchfit(mod_garch, data = returns[,j], solver = "hybrid"),
                     silent=TRUE)
    
    # Save residuals, sigma forecasts, GARCH coefficients, and uniform distribution values
    residuals[, j] <- garch_fit@fit$residuals
    sigma[, j] <- garch_fit@fit$sigma
    garch_coef[[j]] <- garch_fit@fit$coef
    unif_dist[, j] <- fGarch::psstd(q = residuals[, j] / sigma[, j],
                                    nu = garch_coef[[j]][6],
                                    xi = garch_coef[[j]][5])
  }
  
  # Create a list to store the results
  result <- list(residuals = residuals,
                 sigma = sigma,
                 garch_coef = garch_coef,
                 unif_dist = unif_dist)
  
  # Return the results
  return(result)
}



PredictGarch <- function(returns, sigma, zsim, garch_coef) {
  # Initialize matrices to store predicted values
  ret_pred <- matrix(nrow = nrow(zsim), ncol = ncol(returns))
  
  for (j in 1:ncol(returns)) {
    # Get the last observed return and sigma
    ret_t <- as.numeric(returns[nrow(returns), j])
    sigma_t <- sigma[nrow(sigma), j]
    
    # Forecasting at t = 1
    sigma_pred <- sqrt(garch_coef[[j]][7] +  # omega
                         garch_coef[[j]][3] * (ret_t)^2 +  # alpha1
                          garch_coef[[j]][4] * (sigma_t)^2)  # beta1
    mean_pred <- garch_coef[[j]][1] + garch_coef[[j]][2] * ret_t  # mu and ar1
    
    for (i in 1:nrow(zsim)) {
      # building simulated return's matrix
      ret_pred[i, j] <- mean_pred + sigma_pred * zsim[i, j]
    }
  }
  
  return(ret_pred)
}

