LLCG <- function(params, U,
                 copC, copG, copT, copN, copF, copJ){ 
  
  # LLCG: Negative log-likelihood function for estimating copula weights and parameters.
  # Inputs:
  #   params: A numeric vector containing the initial values for copula parameters and weights.
  #   U: A matrix containing the uniform (0 to 1) marginals of the data for each copula.
  #   copC: A copula object (Clayton copula) with initial parameters to be estimated.
  #   copG: A copula object (Gumbel copula) with initial parameters to be estimated.
  #   copT: A copula object (t copula) with initial parameters to be estimated.
  #   copN: A copula object (Normal copula) with initial parameters to be estimated.
  #   copF: A copula object (Frank copula) with initial parameters to be estimated.
  #   copJ: A copula object (Joe copula) with initial parameters to be estimated.
  # Output:
  #   The negative log-likelihood value to be optimized for estimating copula parameters and weights.
  
  
  # Calculate the log-likelihood function to be optimized
  dCop <- matrix(nrow = nrow(U), ncol = 6)
  
  if ("Clayton" %in% names(params)) {
    slot(copC, "parameters") <- params["Clayton"]
    dCop[, 1] <- pi["Clayton"] * copula::dCopula(U, copC)
  } 
  if ("Gumbel" %in% names(params)) {
    slot(copG, "parameters") <- params["Gumbel"]
    dCop[, 2] <- params["piG"] * copula::dCopula(U, copG)
  } 
  if ("t1" %in% names(params)) {
    slot(copT, "parameters") <- c(params["t1"], params["t2"])
    dCop[, 3] <- params["piT"] * copula::dCopula(U, copT)
  } 
  if ("Normal" %in% names(params)) {
    slot(copN, "parameters") <- params["Normal"]
    dCop[, 4] <- params["piN"] * copula::dCopula(U, copN)
  } 
  if ("Frank" %in% names(params)) {
    slot(copF, "parameters") <- params["Frank"]
    dCop[, 5] <- params["piF"] * copula::dCopula(U, copF)
  } 
  if ("Joe" %in% names(params)) {
    slot(copJ, "parameters") <- params["Joe"]
    dCop[, 6] <- params["piJ"] * copula::dCopula(U, copJ)
  } 
  
  # Remove columns with NaN values
  dCop <- dCop[, complete.cases(t(dCop))]
  
  opt <- log(rowSums(dCop))
  
  # Handle infinite values in the log-likelihood
  if(any(is.infinite(opt))){
    opt[which(is.infinite(opt))] <- 0
  }
  
  # Return the negative sum of the log-likelihood
  -sum(opt)
}



eqfun <- function(params, U,
                  copC, copG, copT, copN, copF, copJ){ 
  
  # eqfun: Constrain function to ensure sum of weights = 1.
  # Inputs:
  #   params: A numeric vector containing the initial values for copula parameters and weights.
  #   U: A matrix containing the uniform (0 to 1) marginals of the data for each copula.
  #   copC: A copula object (Clayton copula) with initial parameters to be estimated.
  #   copG: A copula object (Gumbel copula) with initial parameters to be estimated.
  #   copT: A copula object (t copula) with initial parameters to be estimated.
  #   copN: A copula object (Normal copula) with initial parameters to be estimated.
  #   copF: A copula object (Frank copula) with initial parameters to be estimated.
  #   copJ: A copula object (Joe copula) with initial parameters to be estimated.
  # Output:
  #   The sum of the copula weights (pi) to be constrained.
  
  z <- 0
  
  if ("Clayton" %in% names(params)) {
    z <- z + params["piC"]
  } 
  if ("Gumbel" %in% names(params)) {
    z <- z + params["piG"]
  } 
  if ("t1" %in% names(params)) {
    z <- z + params["piT"]
  } 
  if ("Normal" %in% names(params)) {
    z <- z + params["piN"]
  } 
  if ("Frank" %in% names(params)) {
    z <- z + params["piF"]
  } 
  if ("Joe" %in% names(params)) {
    z <- z + params["piJ"]
  } 
  
  return(z)
}



OptMixtureCopulas <- function(unif_dist, K = 10000, combination, pi = NULL) {
  
  # OptMixtureCopulas: Function to optimize the mixture of copulas and generate copula variates.
  # Inputs:
  #   unif_dist: A matrix of uniform marginals (0 to 1) for each copula.
  #   K: Number of copula variates to generate (default = 10000).
  #   combination: Mixture os copulas to optimize - ("Clayton", "Gumbel", "t", "Normal", "Frank", "Joe")
  # Output:
  #   A matrix containing the generated copula variates using the optimized mixture of copulas.
  
  # Initialize copula objects
  ## Creating elliptical copula objects and estimating "initial guesses" for each copula parameter.
  # Then, we maximize log-likelihood of the linear combination of the three copulas
  
  params <- numeric()
  lower <- numeric()
  upper <- numeric()
  
  if(is.null(pi) == TRUE){
    pi = c(Clayton = 1, Gumbel = 1, t = 1, Normal = 1, Frank = 1, Joe = 1)
  }
  
  if ("Clayton" %in% combination) {
    copC <- copula::claytonCopula(param = 10, 
                                  dim = ncol(unif_dist)) # Clayton copula with delta = 2
    params["Clayton"] <- copula::fitCopula(copC, unif_dist, "itau", 
                                           estimate.variance = TRUE)@estimate # Inversion of Kendall's tau for Clayton
    params["piC"] <- pi["Clayton"] 
    # Define lower and upper bounds for the copula parameters and weights
    lower["Clayton"] <- 0.1
    upper["Clayton"] <- copC@param.upbnd
    lower["piC"] <- 0
    upper["piC"] <- 1
  } 
  
  if ("Gumbel" %in% combination) {
    copG <- copula::gumbelCopula(param = 10, 
                                 dim = ncol(unif_dist)) # Gumbel copula with theta = 2
    params["Gumbel"] <- copula::fitCopula(copG, unif_dist, "itau", 
                                          estimate.variance = TRUE)@estimate # Inversion of Kendall's tau for Gumbel
    params["piG"] <- pi["Gumbel"] 
    # Define lower and upper bounds for the copula parameters and weights
    lower["Gumbel"] <- 1
    upper["Gumbel"] <- copG@param.upbnd
    lower["piG"] <- 0
    upper["piG"] <- 1
  } 
  
  if ("t" %in% combination) {
    copT <- copula::tCopula(param = 0.75,
                            dim = ncol(unif_dist),
                            df = ncol(unif_dist))   # t-Copula with parameter 0.5
    params["t1"] <- copula::fitCopula(copT, unif_dist, "mpl", 
                                      estimate.variance = FALSE)@estimate[1] # MPL to estimate Degrees of Freedom (DF)
    params["t2"] <- copula::fitCopula(copT, unif_dist, "mpl", 
                                      estimate.variance = FALSE)@estimate[2] # MPL to estimate Degrees of Freedom (DF)
    params["piT"] <- pi["t"] 
    # Define lower and upper bounds for the copula parameters and weights
    lower["t1"] <- -0.9 
    upper["t1"] <- 1 
    lower["t2"] <- (2 + .Machine$double.eps)
    upper["t2"] <- 100
    lower["piT"] <- 0
    upper["piT"] <- 1
  }
  
  if ("Normal" %in% combination) {
    copN <- copula::normalCopula(param = 0.75, 
                                 dim = ncol(unif_dist)) # Gaussian copula with parameter 0.5
    params["Normal"] <- copula::fitCopula(copN, unif_dist, "mpl", 
                                          estimate.variance = FALSE)@estimate # MPL to estimate parameters for Gaussian copula.
    params["piN"] <- pi["Normal"] 
    # Define lower and upper bounds for the copula parameters and weights
    lower["Normal"] <- copN@param.lowbnd
    upper["Normal"] <- copN@param.upbnd  
    lower["piN"] <- 0
    upper["piN"] <- 1
  } 
  
  if ("Frank" %in% combination) {
    copF <- copula::frankCopula(param = 9, 
                                dim = ncol(unif_dist)) # Frank copula with parameter 1
    params["Frank"] <- copula::fitCopula(copF, unif_dist, "mpl", 
                                         estimate.variance = FALSE)@estimate # MPL to estimate parameters for Frank copula.
    params["piF"] <- pi["Frank"] 
    # Define lower and upper bounds for the copula parameters and weights
    lower["Frank"] <- copF@param.lowbnd
    upper["Frank"] <- copF@param.upbnd 
    lower["piF"] <- 0
    upper["piF"] <- 1
  } 
  
  if ("Joe" %in% combination) {
    copJ <- copula::joeCopula(param = 7,
                              dim = ncol(unif_dist)) # Joe copula
    params["Joe"] <- copula::fitCopula(copJ, unif_dist, "mpl", 
                                       estimate.variance = FALSE)@estimate # MPL to estimate parameters for Joe copula.
    params["piJ"] <- pi["Joe"] 
    # Define lower and upper bounds for the copula parameters and weights
    lower["Joe"] <- copJ@param.lowbnd
    upper["Joe"] <- copJ@param.upbnd  
    lower["piJ"] <- 0
    upper["piJ"] <- 1 
  } 
  
  
  ## Non-linear constrained optimization (RSOLNP)
  opt <- Rsolnp::solnp(pars = params,
                       fun = LLCG,
                       LB = lower,
                       UB = upper,
                       U = unif_dist,
                       eqfun = eqfun,
                       eqB = c(1),
                       copC = copC,
                       copG = copG,
                       copT = copT,
                       copN = copN,
                       copF = copF,
                       copJ = copJ) 
  
  
  ## Generating copula varieties
  # Linear combination of copula varieties
  MixtureCopula <- 0
  if ("Clayton" %in% combination) {
    cC <- opt$pars["piC"] * copula::rCopula(n = K, 
                                            copula = claytonCopula(param = opt$pars["Clayton"],
                                                                   dim = ncol(unif_dist)))
    MixtureCopula <- MixtureCopula + cC
  } 
  
  if ("Gumbel" %in% combination) {
    cG <- opt$pars["piG"] * copula::rCopula(n = K, 
                                            copula = gumbelCopula(param = opt$pars["Gumbel"],
                                                                  dim = ncol(unif_dist)))
    MixtureCopula <- MixtureCopula + cG
  } 
  
  if ("t" %in% combination) {
    cT <- opt$pars["piT"] * copula::rCopula(n = K, 
                                            copula = tCopula(param = opt$pars["t1"],
                                                             df = opt$pars["t2"],
                                                             dim = ncol(unif_dist)))
    MixtureCopula <- MixtureCopula + cT
  }
  
  if ("Normal" %in% combination) {
    cN <- opt$pars["piN"] * copula::rCopula(n = K, 
                                            copula = normalCopula(param = opt$pars["Normal"],
                                                                  dim = ncol(unif_dist)))
    MixtureCopula <- MixtureCopula + cN
  } 
  
  if ("Frank" %in% combination) {
    cF <- opt$pars["piF"] * copula::rCopula(n = K, 
                                            copula = frankCopula(param = opt$pars["Frank"],
                                                                 dim = ncol(unif_dist)))
    MixtureCopula <- MixtureCopula + cF
  } 
  
  if ("Joe" %in% combination) {
    cJ <- opt$pars["piJ"] * copula::rCopula(n = K, 
                                            copula = joeCopula(param = opt$pars["Joe"],
                                                               dim = ncol(unif_dist)))
    MixtureCopula <- MixtureCopula + cJ
  } 
  
  # Return the ctg matrix
  return(MixtureCopula)
}


GaussCopula <- function(unif_dist, K = 10000){
  
  # GaussCopula: Function to generate Gaussian copula variates.
  # Inputs:
  #   unif_dist: A matrix of uniform marginals (0 to 1) for each copula.
  #   K: Number of copula variates to generate (default = 10000).
  # Output:
  #   A matrix containing the generated Gaussian copula variates.
  
  Gcop <- matrix(0, nrow = K, ncol = ncol(unif_dist))
  copn <- copula::normalCopula(param = 0.5, dim = ncol(unif_dist))
  param <- copula::fitCopula(copn, unif_dist, "ml", estimate.variance = FALSE)@estimate 
  Gcop[,] <- copula::rCopula(n = K, 
                             copula = normalCopula(param = param, 
                                                   dim = ncol(unif_dist)))
  return(Gcop)
}



ComputeZSim <- function(copula_mixture, garch_coef) {
  
  # ComputeZSim: Function to compute simulated standardized residuals (zsim).
  # Inputs:
  #   copula_mixture: A matrix containing the copula mixture values.
  #   garch_coef: A list of GARCH coefficients for each copula.
  # Output:
  #   A matrix containing the computed simulated standardized residuals (zsim).
  
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



