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
  opt <- log(pi1 * copula::dCopula(U, copC) + 
               pi2 * copula::dCopula(U, copG) + 
               pi3 * copula::dCopula(U, copt))
  
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


OptMixtureCopulas <- function(unif_dist, K = 10000) {
  # Initialize copula objects
  copt <- copula::tCopula(param = 0.5, dim = ncol(unif_dist))  # t-Copula with parameter 0.5
  copC <- copula::claytonCopula(2, dim = ncol(unif_dist))      # Clayton copula with delta = 2
  copG <- copula::gumbelCopula(2, dim = ncol(unif_dist))       # Gumbel copula with theta = 2
  
  # Define lower and upper bounds for the copula parameters and weights
  lower <- c(0.1, 1, -0.9, (2 + .Machine$double.eps), 0, 0, 0)
  upper <- c(copC@param.upbnd, copG@param.upbnd, 1, 100, 1, 1, 1)
  
  ## Creating elliptical copula objects and estimating "initial guesses" for each copula parameter.
  # Then, we maximize loglikelihood of the linear combination of the three copulas
  par1 <- copula::fitCopula(copC, unif_dist, "itau", estimate.variance = TRUE)@estimate # Inversion of Kendall's tau for Clayton
  par2 <- copula::fitCopula(copG, unif_dist, "itau", estimate.variance = TRUE)@estimate # Inversion of Kendall's tau for Gumbel
  par3 <- copula::fitCopula(copt, unif_dist, "mpl", estimate.variance = FALSE)@estimate # MPL to estimate Degrees of Freedom (DF)
  
  # Initialize weights for copulas (initial guesses = 1/3 each)
  par6 <- par5 <- par4 <- 1/3
  
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
  ctg <- Cc <- Cg <- Ct <- matrix(nrow = K, ncol = ncol(unif_dist))
  
  ## Generating copula variates
  Cc[, ] <- cop_param[5] * copula::rCopula(n = K,
                                           copula = copula::claytonCopula(param = cop_param[1],
                                                                          dim = ncol(unif_dist)))
  Cg[, ] <- cop_param[6] * copula::rCopula(n = K,
                                           copula = copula::gumbelCopula(param = cop_param[2],
                                                                         dim = ncol(unif_dist)))
  Ct[, ] <- cop_param[7] * copula::rCopula(n = K,
                                           copula = copula::tCopula(param = cop_param[3],
                                                                    df = cop_param[4],
                                                                    dim = ncol(unif_dist)))
  
  # Linear combination of copula varieties
  ctg <- Cc + Ct + Cg
  
  # Return the ctg matrix
  return(ctg)
}


GaussCopula <- function(unif_dist, K = 10000){
  Gcop <- matrix(0, nrow = K, ncol = ncol(unif_dist))
  Gcop[,] <- copula::rCopula(n = K, 
                             copula = normalCopula(param = 0.5, 
                                                   dim = ncol(unif_dist)))
  return(Gcop)
}


ComputeZSim <- function(copula_mixture, garch_coef) {
  # Create an empty matrix to store zsim values
  zsim <- matrix(nrow = nrow(copula_mixture), ncol = ncol(copula_mixture))
  
  # Compute zsim values for each column of copula_mixture
  for (j in 1:ncol(copula_mixture)) {
    # Compute zsim using qsstd function from fGarch package
    zsim[, j] <- fGarch::qsstd(copula_mixture[, j],
                               nu = garch_coef[[j]][6],
                               xi = garch_coef[[j]][5]) /
      sd(fGarch::qsstd(copula_mixture[, j], nu = garch_coef[[j]][6],
                       xi = garch_coef[[j]][5]))
  }
  
  # Return the zsim matrix
  return(zsim)
}



