#################################################################################
# Main file to run worst-case CVaR portfolio optimization based mixture copulas #
# Authors: Joao Ramos Jungblut ##################################################
# Last update: 2023-07-24 #######################################################
#################################################################################


# setting R project environment
my_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(my_dir)


# cleaning variables and graphs
rm(list=ls())
graphics.off()


# Load required packages
library(tidyverse)     # Data manipulation and visualization
library(tidyquant)     # Financial data analysis
library(readxl)        # Read excel files
library(rugarch)       # Univariate GARCH modeling
library(fGarch)        # Multivariate GARCH modeling
library(copula)        # Copula modeling
library(Rsolnp)        # Nonlinear optimization
library(fPortfolio)    # Portfolio optimization
library(PerformanceAnalytics) # Performance metrics
library(xts)           # Time series object
library(timetk)        # Time series object
library(xtable)        # Create LaTex tables
library(ggplot2)       # Produce graph
library("ROI")         # Optimization
library("ROI.plugin.glpk") # Optimization
library("ROI.plugin.quadprog") # Optimization
library("ROI.plugin.alabama") # Optimization


# Fetch data
# OBSERVATION: CSIP6 was removed from data due to impossibility of found data
Ret <- read_excel("data_directory/StockPrice.xlsx")[-1,] %>% 
  mutate(date = as.Date(data), # convert to date format
         across(VALE3:WHMT3, as.numeric)) %>%  # convert to numeric format
  select(date, everything(), -data) %>% # select date and tickers
  gather("Ticker", "AdjClose", -date, na.rm = TRUE) %>% # adjust to panel data
  group_by(Ticker) %>% # calculate return for each ticker
  mutate(Ret = log(AdjClose/ dplyr::lag(AdjClose))) %>%  # compute log return
  ungroup() %>% # stop grouping tickers
  select(date, Ticker, Ret) %>% # selecting necessary columns
  spread(key = "Ticker", value = "Ret", fill = 0) # transforming data 
Ret


# Choosing stocks on Ibovespa per year
MarketComposition <- read_excel("data_directory/MarketComposition.xlsx") %>% 
  mutate(date = as.Date(Data), # convert to date format
         m = month(date, label = TRUE)) %>% # convert to numeric format
  dplyr::filter(m == "dez") %>% # update factors 
  select(date, everything(),-c(Data, m)) %>% # select date, month and tickers
  gather("Ticker", "Composition", -date, na.rm=TRUE) %>% # adjust to panel data
  group_by(date) %>% # Select tickers by last date of the month
  reframe(Tickers = Ticker)
MarketComposition

Update <- MarketComposition %>% 
  select(date) %>% 
  unique() %>% # avoid repetition
  deframe() # create vector
Update

Symbols <- map(Update, .f = function(x){ # list of vectors with Ibov tickers per year
  MarketComposition %>% 
    dplyr::filter(date == x) %>% # choose stocks presented on Ibov
    select(Tickers) %>% 
    unique() %>% # avoid repetition
    deframe() # create vector
})
names(Symbols) <- Update
Symbols


# Filtering stocks on Ibovespa and splitting  data set
Ret_inSample <- map(Update, .f = function(x){ 
  Ret %>% 
    select(date,  Symbols[[paste(x)]]) %>% # filter Ibov stocks
    dplyr::filter(date > as.Date(x) - 365, # split data in in-sample
                  date <= as.Date(x)) %>% 
    select(-date) %>% # exclude date column
    as.matrix() # convert to matrix format
})
names(Ret_inSample) <- Update
Ret_inSample

Ret_outofSample <- map(Update, .f = function(x){ 
  Ret %>% 
    select(date,  Symbols[[paste(x)]]) %>% # filter Ibov stocks
    dplyr::filter(date > as.Date(x), # split data in out-of-sample
                  date <= as.Date(x) + 365) %>% 
    select(-date) %>% # exclude date column
    as.matrix() # convert to matrix format
})
names(Ret_outofSample) <- Update + 365
Ret_outofSample


# Initialize the GARCH specification
mod_garch <- try(rugarch::ugarchspec(variance.model = list(model = "sGARCH", 
                                                           garchOrder = c(1, 1),
                                                           variance.targeting = TRUE),
                                     mean.model = list(armaOrder = c(1, 0)),
                                     distribution.model = "sstd"),
                 silent = TRUE)

# Create returns matrix
returns <- Ret_inSample[["1997-12-31"]]


# Matrix to save model coefficients, standard deviation and residuals
unif_dist <- garch_pred <- sigma <- residuals <- matrix(nrow = nrow(returns), 
                                                        ncol = ncol(returns))
garch_coef <- vector("list", length = ncol(returns))
colnames(unif_dist) <- colnames(returns)
colnames(garch_pred) <- colnames(returns)
colnames(sigma) <- colnames(returns)
colnames(residuals) <- colnames(returns)
names(garch_coef) <- colnames(returns)


# Loop over each asset
for(j in 1:ncol(returns)){
  # Fitting GARCH model
  print(j)
  garch_fit <- try(ugarchfit(mod_garch, data = returns[,j], solver = "hybrid"),
                   silent=TRUE)
  
  # Check if garch_fit is an error or a successful fit
  if (inherits(garch_fit, "try-error")) {
    # Handle the error (e.g., print an error message)
    cat("Error in fitting GARCH model for asset", j, "\n")
  } else {
    # Save residuals, sigma, GARCH coefficients, and uniform distribution values
    residuals[, j] <- garch_fit@fit$residuals
    sigma[, j] <- garch_fit@fit$sigma
    garch_coef[[j]] <- garch_fit@fit$coef
    unif_dist[, j] <- fGarch::psstd(q = residuals[, j] / sigma[, j],
                                    nu = garch_coef[[j]][6],
                                    xi = garch_coef[[j]][5])
  }
}


# Identify columns with complete cases 
complete_columns <- apply(unif_dist, 2, function(x) all(!is.na(x)))


# Subset the matrix to keep only columns with complete cases
unif_dist <- unif_dist[, complete_columns] # drop na columns


# Initialize copula objects
copT <- copula::tCopula(param = 0.5,
                        dim = ncol(unif_dist),
                        df = ncol(unif_dist))   # t-Copula with parameter 0.5
copC <- copula::claytonCopula(param = 2, 
                              dim = ncol(unif_dist)) # Clayton copula with delta = 2
copG <- copula::gumbelCopula(param = 2, 
                             dim = ncol(unif_dist)) # Gumbel copula with theta = 2
copF <- copula::frankCopula(param = 1, 
                            dim = ncol(unif_dist)) # Frank copula with parameter 1
copJ <- copula::joeCopula(param = 2,
                          dim = ncol(unif_dist)) # Joe copula
copN <- copula::normalCopula(param = 0.5, 
                             dim = ncol(unif_dist)) # Gaussian copula with parameter 0.5


# Define lower and upper bounds for the copula parameters and weights
lower <- c(0.1, 1, -0.9, (2 + .Machine$double.eps), 0, 0, 0)
upper <- c(copC@param.upbnd, copG@param.upbnd, 1, 100, 1, 1, 1)


## Creating elliptical copula objects and estimating "initial guesses" for each copula parameter.
# Then, we maximize log-likelihood of the linear combination of the three copulas
paramC <- copula::fitCopula(copC, unif_dist, "itau", estimate.variance = TRUE)@estimate # Inversion of Kendall's tau for Clayton
paramG <- copula::fitCopula(copG, unif_dist, "itau", estimate.variance = TRUE)@estimate # Inversion of Kendall's tau for Gumbel
paramT <- copula::fitCopula(copT, unif_dist, "mpl", estimate.variance = FALSE)@estimate # MPL to estimate Degrees of Freedom (DF)
paramN <- copula::fitCopula(copN, unif_dist, "mpl", estimate.variance = FALSE)@estimate
paramF <- copula::fitCopula(copF, unif_dist, "mpl", estimate.variance = FALSE)@estimate
paramJ <- copula::fitCopula(copJ, unif_dist, "mpl", estimate.variance = FALSE)@estimate


# Create params vector
params <- c(Clayton = paramC, 
            Gumbel = paramG,
            t = paramT,
            Gaussian = paramN,
            Frank = paramF,
            Joe = paramJ)


LLCG <- function(params, U){ 
  
  # LLCG: Negative log-likelihood function for estimating copula weights and parameters.
  # Inputs:
  #   params: A numeric vector containing the initial values for copula parameters and weights.
  #   U: A matrix containing the uniform (0 to 1) marginals of the data for each copula.
  # Output:
  #   The negative log-likelihood value to be optimized for estimating copula parameters and weights.
  
  
  
  if ("Clayton" %in% names(params)) {
    slot(copC, "parameters") <- params["Clayton"]
  } 
  if ("Gumbel" %in% names(params)) {
    slot(copG, "parameters") <- params["Gumbel"]
  } 
  if ("t" %in% names(params)) {
    slot(copT, "parameters") <- params["t"]
  } 
  if ("Gaussian" %in% names(params)) {
    slot(copN, "parameters") <- params["Gaussian"]
  } 
  if ("Frank" %in% names(params)) {
    slot(copF, "parameters") <- params["Frank"]
  } 
  if ("Joe" %in% names(params)) {
    slot(copJ, "parameters") <- params["Joe"]
  } 
  
  
  # Set copula weights
  pi <- 1/length(params)
  
  
  # Calculate the log-likelihood function to be optimized
  dCop <- matrix(nrow = nrow(U), ncol = length(params))
  colnames(dCop) <- names(params)
  
  for (i in seq_along(params)) {
    if ("Clayton" %in% names(params)) {
      dCop[, i] <- pi * copula::dCopula(U, copC)
    } else if ("Gumbel" %in% names(params)) {
      dCop[, i] <- pi * copula::dCopula(U, copG)
    } else if ("t" %in% names(params)) {
      dCop[, i] <- pi * copula::dCopula(U, copT)
    } else if ("Gaussian" %in% names(params)) {
      dCop[, i] <- pi * copula::dCopula(U, copN)
    } else if ("Frank" %in% names(params)) {
      dCop[, i] <- pi * copula::dCopula(U, copF)
    } else if ("Joe" %in% names(params)) {
      dCop[, i] <- pi * copula::dCopula(U, copJ)
    }
  }
  opt <- log(rowSums(dCop))
  
  # Handle infinite values in the log-likelihood
  if(any(is.infinite(opt))){
    opt[which(is.infinite(opt))] <- 0
  }
  
  # Return the negative sum of the log-likelihood
  -sum(opt)
}


# Initialize lower and upper bounds as named vectors
lower <- numeric()
upper <- numeric()

if ("Clayton" %in% names(params)) {
  lower["Clayton"] <- 0.1
  upper["Clayton"] <- copC@param.upbnd
} 
if ("Gumbel" %in% names(params)) {
  lower["Gumbel"] <- 1
  upper["Gumbel"] <- copG@param.upbnd
}
if ("t1" %in% names(params)) {
  lower["t1"] <- -0.9 
  upper["t1"] <- 1 
} 
if ("t2" %in% names(params)) {
  lower["t2"] <- (2 + .Machine$double.eps)
  upper["t2"] <- 100
} 
if ("Gaussian" %in% names(params)) {
  lower["Gaussian"] <- copN@param.lowbnd
  upper["Gaussian"] <- copN@param.upbnd  
} 
if ("Frank" %in% names(params)) {
  lower["Frank"] <- copF@param.lowbnd
  upper["Frank"] <- copF@param.upbnd  
} 
if ("Joe" %in% names(params)) {
  lower["Joe"] <- copJ@param.lowbnd
  upper["Joe"] <- copJ@param.upbnd  
}


## Non-linear constrained optimization (RSOLNP)
opt <- Rsolnp::solnp(pars = params,
                     fun = LLCG,
                     LB = lower,
                     UB = upper,
                     U = unif_dist)


## Generating copula variaties
cC <- copula::rCopula(n = 10000, copula = copC)
cG <- copula::rCopula(n = 10000, copula = copG)
ct <- copula::rCopula(n = 10000, copula = copt)
cn <- copula::rCopula(n = 10000, copula = copn)
cF <- copula::rCopula(n = 10000, copula = copF)
cJ <- copula::rCopula(n = 10000, copula = copJ)

plot(ct)



