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


# Importing modules
source("data_preprocessing.R")
source("garch_estimate.R")
source("copula_estimate.R")
source("portfolio_optimization.R")
source("portfolio_analysis.R")
source("performance_metrics.R")
source("exporting_results.R")


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


# Define a function to perform computations for each  window
Pipeline <- function(inSample, outofSample, Update, copulas){
  
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
                                Alpha = 0.05, 
                                TargetReturn = 0,
                                NumAssets = 16)
    
    # Calculate portfolio returns based on the optimal weights 
    ret_matrix_outofsample <- outofSample[[paste(as.Date(x) + 365)]][,colnames(returns)] # select valid stocks
    portfolio_returns <- ret_matrix_outofsample  %*%  weights
  }) %>% 
    bind_rows()
  
  return(Pipe)
}


# Calculate cumulative returns
results <- Pipeline(Ret_inSample, Ret_outofSample, Update, copulas = c("Clayton", "Joe"))
cumulative_returns <- cumprod(1 + results) - 1


# Plot the cumulative returns
plot(cumulative_returns, type = "l", col = "green", lwd = 1,
     main = "Portfolio Performance",
     xlab = "Time Period", ylab = "Cumulative Returns")

















