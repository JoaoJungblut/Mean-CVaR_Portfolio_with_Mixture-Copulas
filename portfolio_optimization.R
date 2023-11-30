

CVaROptimization <- function(returns, Alpha = 0.025, TargetReturn = 0.0003) { # Target return must consider transaction costs
  
  # CVaROptimization: Function to perform CVaR optimization for portfolio weights.
  # Inputs:
  #   returns: A matrix of asset returns.
  #   Alpha: CVaR alpha level (default = 0.025).
  #   TargetReturn: Daily target return constraint (default = 0.0003).
  # Output:
  #   A matrix containing the optimized portfolio weights using CVaR optimization.
  
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



CVaROptimization <- function(returns, 
                             Alpha = 0.975, 
                             TargetReturn = 0,
                             Turnover = 0.0003,
                             NumAssets = 8) {
  
  # CVaROptimization: Function to perform CVaR optimization for portfolio weights.
  # Inputs:
  #   returns: A data frame of asset returns.
  #   Alpha: CVaR alpha level (default = 0.025).
  #   TargetReturn: Daily target return constraint (default = 0).
  #   Turnover: Daily turnover measure constraint (default = 0.0003) equal to transaction costs. 
  #   NumAssets: Maximum number of assets in the portfolio (default = 8).
  # Output:
  #   A matrix containing the optimized portfolio weights using CVaR optimization.
  
  cvar_objective <- function(r_mat, alpha, probs = NULL) {
    x.names <- colnames(r_mat)
    N <- NCOL(r_mat)
    S <- NROW(r_mat)
    mu <- colMeans(r_mat)
    if (is.null(probs)) probs <- rep(1/S, S)
    if (alpha < 0.5) alpha <- 1 - alpha
    
    Amat <- cbind(as.matrix(r_mat),  diag(S), 1)
    var.names <- c(x.names, paste0("z_cvar_aux", seq_len(S)), "gamma")
    
    ## set bounds for gamma (-Inf, Inf) 
    bnds <- ROI::V_bound(li = c(N + S + 1), lb = c( -Inf),
                         ui = c(N + S + 1), ub = c(  Inf))
    
    constraint <- L_constraint(L = Amat, dir = rep(">=", S), 
                               rhs = rep(0, S), 
                               names = var.names)
    
    objective <- L_objective(c(rep(0, N), probs/(1 - alpha), 1))
    
    list(objective = objective, constraint = constraint, bounds = bnds)
  }
  
  budget_constraint <- function(r_mat, dir = "==", rhs = 1) {
    x.names <- colnames(r_mat)
    L_constraint(L = rep(1, NCOL(r_mat)), 
                 dir = dir,  rhs = rhs, names = x.names)
  }
  
  reward_constraint <- function(r_mat, dir = ">=", rhs = TargetReturn) {
    x.names <- colnames(r_mat)
    L_constraint(L = colMeans(r_mat), dir = dir,  
                 rhs = rhs, names = x.names)
  }
  
  turnover_constraint <- function(r_mat, x0 = NULL, dir = "<=", rhs = Turnover) {
    x.names <- colnames(r_mat)
    N <- NCOL(r_mat)
    S <- NROW(r_mat)
    if (is.null(x0)) x0 <- rep(1/N, N)
    Amat <- cbind(diag(N), - diag(N), diag(N))
    var.names <- c(x.names,  
                   paste0("y_plus_aux", seq_len(N)), 
                   paste0("y_minus_aux", seq_len(N)))
    
    rbind(L_constraint(L = Amat, dir = rep("==", N), rhs = x0, 
                       names = var.names),
          L_constraint(c(rep(0, N), rep(1, N), rep(1, N) ), 
                       dir = dir, rhs = rhs, names = var.names))
  }
  
  cardinality_constraint <- function(r_mat, dir = "<=", rhs = NumAssets) {
    x.names <- colnames(r_mat)
    N <- NCOL(r_mat)
    Amat <- cbind(diag(N), -diag(N))
    var.names <- c(x.names, paste0("z_card_aux", seq_len(N)))
    cat("Variable types for z_card_aux must be set to binary.\n")
    rbind(L_constraint(L = Amat, dir = rep("<=", N),
                       rhs = rep(0, N), names = var.names),
          L_constraint(L = c(rep(0, N), rep(1, N)), dir = dir,
                       rhs = rhs, names = var.names))
  }
  
  tmp <- cvar_objective(returns, 0.975)
  lp <- OP()
  
  constraints(lp) <- rbind(
    tmp$constraint,
    budget_constraint(returns),
    reward_constraint(returns),
    #turnover_constraint(returns), 
    cardinality_constraint(returns, dir = "<=", rhs = NumAssets),
    use.names = TRUE
  )
  
  ## Variable types for z_card_aux must be set to binary.
  obj <- c((tmp$objective)$L)
  objective(lp) <- c(obj, double(NCOL(constraints(lp)) - length(obj)))
  ## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'constraints' for signature '"OP"'
  types(lp) <- rep("C",  NCOL(constraints(lp)))
  ## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'constraints' for signature '"OP"'
  types(lp)[grep("z_card_aux", constraints(lp)$names)] <- "B"
  ## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'constraints' for signature '"OP"'
  (sol <- ROI_solve(lp, solver = "glpk"))
  ## Error in ROI_solve(lp, solver = "glpk"): objective is missing, with no default
  weights <- round(solution(sol)[1:ncol(returns)], 4)[1:ncol(returns)]
  
  return(weights)
}


CVaROptimization <- function(returns, 
                             Alpha = 0.95, 
                             TargetReturn = 0,
                             NumAssets = 8) {
  
  # CVaROptimization: Function to perform CVaR optimization for portfolio weights.
  # Inputs:
  #   returns: A data frame of asset returns.
  #   Alpha: CVaR alpha level (default = 0.025).
  #   TargetReturn: Daily target return constraint (default = 0).
  #   Turnover: Daily turnover measure constraint (default = 0.0003) equal to transaction costs. 
  #   NumAssets: Maximum number of assets in the portfolio (default = 8).
  # Output:
  #   A matrix containing the optimized portfolio weights using CVaR optimization.
  
  m <- model()
  m$variable(portfolio, lb = 0)
  m$minimize(cvar(portfolio, alpha = Alpha))
  m$subject_to(budget_norm(portfolio))
  m$subject_to(reward(portfolio) >= TargetReturn)
  m$subject_to(cardinality(portfolio) <= NumAssets)
  opt <- optimize(m, solver="glpk", data=list(returns = as.matrix(returns)))
  weights <- round(opt$solution[grep("portfolio", names(opt$solution))], 3)
  
  return(weights)
}



NaiveDiversification <- function(returns) {
  
  # NaiveDiversification: Function to calculate portfolio returns using naive diversification (equal weights).
  # Inputs:
  #   returns: A data frame containing asset returns with a 'date' column and individual asset columns.
  # Output:
  #   A data frame with 'date' and 'portfolio_return' columns representing portfolio returns.
  
  # Calculate the number of assets
  num_assets <- ncol(returns) - 1  # Subtract 1 for the 'date' column
  
  # Calculate the equal weights for each asset
  weight_per_asset <- 1 / num_assets
  
  # Extract the returns data without the 'date' column
  returns_df <- returns[, -1]  # Drop the 'date' column
  
  # Calculate portfolio returns for each date
  portfolio_returns <- (rowSums(returns_df, na.rm = TRUE) * weight_per_asset) - 0.0003 # minus transaction costs
  
  # Create a new data frame with the date and portfolio returns
  portfolio_returns_df <- data.frame(date = returns$date, portfolio_return = portfolio_returns)
  
  return(portfolio_returns_df)
}

