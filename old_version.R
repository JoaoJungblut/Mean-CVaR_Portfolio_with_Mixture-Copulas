#################################################################################
# Main file to run worst-case CVaR portfolio optimization based mixture copulas #
# Authors: Joao Ramos Jungblut ##################################################
# Last update: 2023-06-27 #######################################################
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
library(rugarch)       # Univariate GARCH modeling
library(fGarch)        # Multivariate GARCH modeling
library(copula)        # Copula modeling
library(Rsolnp)        # Nonlinear optimization
library(fPortfolio)    # Portfolio optimization
library(PerformanceAnalytics) # Performance metrics
library(xts) # Time series object

# Importing modules
source("data_preprocessing.R")
source("garch_estimate.R")
source("copula_estimate.R")
source("portfolio_optimization.R")
source("performance_metrics.R")

# Define the list of stock tickers and the start date for data retrieval
#tickers <- c("PETR4.SA", "VALE3.SA", "ITUB4.SA", "BBAS3.SA", "ABEV3.SA", 
#             "BBDC4.SA", "GRND3.SA", "SMTO3.SA", "SLCE3.SA", "VIVT3.SA")
#start_date <- "2010-01-01"

# Retrieve the stock returns for the given tickers and start date
#returns <- GetReturns(tickers = tickers, start_date = start_date)
returns <- read_csv("data_directory/etfs_rtn.csv")[-1]
###### Error in estimating parameters of Copula t for > 28 assets

# Creating auxiliary matrices and list
N <- base::ncol(returns) - 1   # Number of assets
We <- 252 # Window size for GARCH estimation
Wt <- nrow(returns) # Total size of window 
K <- 1000 # Number of simulations 
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

set.seed(64)

# Rolling window estimation
for (i in (We + 1):Wt){
  print(paste(i - We, "of", Wt - We))
  
  # Establishing window interval in-sample
  t1 <- i - We
  t2 <- i - 1
  
  # Convert the in-sample returns data to a matrix format
  ret_matrix_insample <- as.matrix(returns[t1:t2, -1])
  
  # Create a logical vector indicating if each asset has sufficient data
  assets_with_valid_returns <- !colMeans(is.na(ret_matrix_insample[,]))
  
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
  portfolio_returns[i,] <- RetPortfolio(returns = ret_matrix_outofsample,  
                                        weights = rbind(weights[i,])) - 0.0003 # minus the transaction costs
  
  benchmark_portfolio_returns[i,] <- RetPortfolio(returns = ret_matrix_outofsample,  
                                                  weights = rbind(weights_benchmark[i,])) - 0.0003 # minus the transaction cost

}

# Construct naive diversification portfolio 
naive_portfolio <- NaiveDiversification(returns)

# Convert the portfolio_returns matrix to an xts object
portfolio_returns_xts <- xts::xts(portfolio_returns[We:Wt,], 
                                  order.by = returns[We:Wt,]$date)

benchmark_returns_xts <- xts::xts(benchmark_portfolio_returns[We:Wt,],
                                  order.by = returns[We:Wt,]$date)

naive_returns_xts <- xts::xts(naive_portfolio[We:Wt, -1],
                              order.by = naive_portfolio[We:Wt,]$date)

# Compute Performance
portfolio_performance <- ComputePerformance(portfolio_returns_xts)
benchmark_performance <- ComputePerformance(benchmark_returns_xts)
naive_performance <- ComputePerformance(naive_returns_xts)

# Generate graph
PerformanceAnalytics::charts.PerformanceSummary(portfolio_returns_xts)
PerformanceAnalytics::charts.PerformanceSummary(benchmark_returns_xts)
PerformanceAnalytics::charts.PerformanceSummary(naive_returns_xts)
