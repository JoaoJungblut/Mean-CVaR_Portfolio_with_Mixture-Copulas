# Load required packages
library(tidyverse)     # Data manipulation and visualization
library(tidyquant)     # Financial data analysis
library(rugarch)       # Univariate GARCH modeling
library(fGarch)        # Multivariate GARCH modeling
library(copula)        # Copula modeling
library(Rsolnp)        # Nonlinear optimization
library(fPortfolio)    # Portfolio optimization

# Importing data and calculating returns
returns <- tidyquant::tq_get(c("PETR4.SA", "VALE3.SA", "ITUB4.SA",
                               "BBAS3.SA", "ABEV3.SA", "BBDC4.SA",
                               "SLCE3.SA", "BBSE3.SA", "SMTO3.SA",
                               "AMER3.SA", "MGLU3.SA", "KLBN11.SA",
                               "ITSA4.SA", "BPAC11.SA", "SUZB3.SA",
                               "RENT3.SA", "LREN3.SA", "PETR3.SA",
                               "SANB11.SA", "B3SA3.SA", "YDUQ3.SA"), 
                             from = "1997-01-01") %>% 
  dplyr::select(date, symbol, adjusted) %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::mutate(return = log(adjusted) - log(dplyr::lag(adjusted))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-adjusted) %>% 
  tidyr::spread(symbol, return) %>% 
  na.omit()

# Creating auxiliary matrices and list
N <- base::ncol(returns) - 1   # Number of assets
K <- 252                 # Window size for GARCH estimation
index_vector <- seq(1, nrow(returns), by = K)  # Index vector for rolling optimization
names_vector <- names(returns)[-1]   # Asset names for reference
portfolio_daily_returns <- matrix(nrow = nrow(returns), ncol = 1)  # Matrix to store portfolio returns
portfolio_daily_returns[1:K, ] <- 0  # Initialize the first K rows as zero

# Initialize the GARCH specification
mod_garch <- try(
  rugarch::ugarchspec(variance.model = list(model = "sGARCH",
                                            garchOrder = c(1, 1),
                                            variance.targeting = TRUE),
                      mean.model = list(armaOrder = c(1, 0)),
                      distribution.model = "sstd"),
  silent = TRUE
)

# Initialize copula objects
copt <- copula::tCopula(param = 0.5, dim = N)  # t-Copula with parameter 0.5
copC <- copula::claytonCopula(2, dim = N)      # Clayton copula with delta = 2
copG <- copula::gumbelCopula(2, dim = N)       # Gumbel copula with theta = 2

# Define lower and upper bounds for the copula parameters and weights
lower <- c(0.1, 1, -0.9, (2 + .Machine$double.eps), 0, 0, 0)     
upper <- c(copC@param.upbnd, copG@param.upbnd, 1, 100, 1, 1, 1) 

# Initialize weights for copulas (initial guesses = 1/3 each)
par6 <- par5 <- par4 <- 1/3 

# Negative log-likelihood function for estimating copula weights and parameters
LLCG <- function(params, U, copC, copG, copt){ 
  # Set copula parameters
  slot(copC, "parameters") <- params[1]    # Initial Clayton parameter
  slot(copG, "parameters") <- params[2]    # Initial Gumbel parameter 
  slot(copt, "parameters") <- params[3:4]  # Initial t parameters (correlation and degrees of freedom)
  
  # Set copula weights
  pi1 <- params[5]  # Weight of Clayton copula
  pi2 <- params[6]  # Weight of Gumbel copula
  pi3 <- params[7]  # Weight of t copula
  
  # Calculate the log-likelihood function to be optimized
  opt <- log(pi1 * copula::dCopula(U, copC) + pi2 * copula::dCopula(U, copG) + pi3 * copula::dCopula(U, copt))
  
  # Handle infinite values in the log-likelihood
  if(any(is.infinite(opt))){
    opt[which(is.infinite(opt))] <- 0
  }
  
  # Return the negative sum of the log-likelihood
  -sum(opt)
}

# Constrain function to ensure sum(weights) = 1
eqfun <- function(params, U, copC, copG, copt){ 
  z <- params[5] + params[6] + params[7]
  return(z)
}

##### Setting up a fPortfolio min CVaR optimization
cvar_opt <- matrix(nrow = length(index_vector) - 1, ncol = N)  # Matrix to store CVaR optimization results
targetReturn <- 0  # Daily target return constraint
frontierSpec <- fPortfolio::portfolioSpec()  # Portfolio specification for optimization
fPortfolio::setType(frontierSpec) <- "CVaR"  # Set portfolio type as CVaR
fPortfolio::setSolver(frontierSpec) <- "solveRglpk.CVAR"  # Use linear programming solver for CVaR optimization
fPortfolio::setAlpha(frontierSpec) <- 0.025  # Set CVaR alpha level as 0.05 (CVaR_0.95)
fPortfolio::setTargetReturn(frontierSpec) <- targetReturn  # Set the daily target return constraint

# Iterating over n assets
for (i in 2:length(index_vector)) {
  
  # Matrix to save sigma forecasts
  unif_dist <- garch_pred <- sigma <- residuals <- matrix(nrow = K, ncol = N)
  garch_coef <- vector("list", length = N)
  
  t1 <- index_vector[i - 1]
  t2 <- index_vector[i] - 1
  
  # Looping through each asset
  for (j in 1:length(names_vector)) {
    x <- cbind(returns[t1:t2, (j + 1)])
    
    # Fitting GARCH model
    garch_fit <- try(rugarch::ugarchfit(mod_garch, data = x, solver = "hybrid"),
                     silent = TRUE)
    
    residuals[, j] <- garch_fit@fit$residuals
    sigma[, j] <- garch_fit@fit$sigma
    garch_coef[[j]] <- garch_fit@fit$coef
    unif_dist[, j] <- fGarch::psstd(q = residuals[, j] / sigma[, j],
                                    nu = garch_coef[[j]][6],
                                    xi = garch_coef[[j]][5])
  }
  
  ## Creating elliptical copula objects and estimating "initial guesses" for each copula parameter.
  # Then, we maximize loglikelihood of the linear combination of the three copulas
  par1 <- copula::fitCopula(copC, unif_dist, "itau", estimate.variance = TRUE)@estimate # Inversion of Kendall's tau for Clayton
  par2 <- copula::fitCopula(copG, unif_dist, "itau", estimate.variance = TRUE)@estimate # Inversion of Kendall's tau for Gumbel
  par3 <- copula::fitCopula(copt, unif_dist, "mpl", estimate.variance = FALSE)@estimate # MPL to estimate Degrees of Freedom (DF)
  
  ## Non-linear constrained optimization (RSOLNP)
  opt <- Rsolnp::solnp(pars = c(par1, par2, par3, par4, par5, par6),
                       fun = LLCG,
                       LB = lower,
                       UB = upper,
                       copt = copt,
                       copC = copC,
                       copG = copG,
                       U = unif_dist,
                       eqfun = eqfun,
                       eqB = c(1)) # RSOLNP
  
  ## Saving optimization parameters in a list
  cop_param <- opt$pars
  
  # Clayton, t, gumbel, and ctg variates matrix
  ctg <- Cc <- Cg <- Ct <- matrix(nrow = K, ncol = N)
  
  ## Generating copula variates
  Cc[, ] <- cop_param[5] * copula::rCopula(n = K,
                                           copula = copula::claytonCopula(param = cop_param[1],
                                                                          dim = N))
  Cg[, ] <- cop_param[6] * copula::rCopula(n = K,
                                           copula = copula::gumbelCopula(param = cop_param[2],
                                                                         dim = N))
  Ct[, ] <- cop_param[7] * copula::rCopula(n = K,
                                           copula = copula::tCopula(param = cop_param[3],
                                                                    df = cop_param[4],
                                                                    dim = N))
  
  # Linear combination of copula varieties
  ctg <- Cc + Ct + Cg
  
  # For each asset, generate copula 'z' dependence structure
  rtn_pred <- mean_pred <- sigma_pred <- zsim <- matrix(nrow = K, ncol = N)
  for (j in 1:N) {
    
    zsim[, j] <- fGarch::qsstd(ctg[, j],
                               nu = garch_coef[[j]][6],
                               xi = garch_coef[[j]][5]) /
      sd(fGarch::qsstd(ctg[, j], nu = garch_coef[[j]][6],
                       xi = garch_coef[[j]][5]))
    
    rtn_t <- as.numeric(returns[(t2 - 1), (j + 1)])
    sigma_t <- sigma[K, j]
    
    sigma_pred[1, j] <- sqrt(garch_coef[[j]][7] +  # omega
                               garch_coef[[j]][3] * (rtn_t)^2 +  # alpha1
                               garch_coef[[j]][4] * (sigma_t)^2)  # beta1
    mean_pred[1, j] <- garch_coef[[j]][1] + garch_coef[[j]][2] * rtn_t  # mu and ar1
    rtn_pred[1, j] <- mean_pred[1, j] + sigma_pred[1, j] * zsim[1, j]
    
    for (t in 2:K) {
      sigma_pred[t, j] <- sqrt(garch_coef[[j]][7] +  # omega
                                 garch_coef[[j]][3] * (rtn_pred[(t - 1), j])^2 +  # alpha1
                                 garch_coef[[j]][4] * (sigma_pred[(t - 1), j])^2)  # beta1
      mean_pred[t, j] <- garch_coef[[j]][1] + garch_coef[[j]][2] * rtn_pred[(t - 1), j]  # mu and ar1
      rtn_pred[t, j] <- mean_pred[t, j] + sigma_pred[t, j] * zsim[t, j]
    }
  }
  
  ## Optimizing portfolio using K simulated returns for each asset, for optimization period i
  returnfPort <- as.timeSeries(rtn_pred[, 1:N])
  frontier1g <- fPortfolio::efficientPortfolio(data = returnfPort,
                                               spec = frontierSpec,
                                               constraints = "LongOnly")
  cvar_opt[(i - 1), 1:N] <- fPortfolio::getWeights(frontier1g)  # Storing resulting weights
  
  # Calculate portfolio returns
  t1 <- t1+K
  t2 <- min(nrow(returns), t2+K)
  portfolio_returns <- as.matrix(returns[t1:t2, names_vector]) %*% t(cvar_opt[(-1), ])
  portfolio_daily_returns[t1:t2, ] <- rowSums(portfolio_returns)
  
}


