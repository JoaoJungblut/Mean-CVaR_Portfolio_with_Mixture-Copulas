#################################################################################
# Main file to run worst-case CVaR portfolio optimization based mixture copulas #
# Authors: Joao Ramos Jungblut ##################################################
# Last update: 2023-07-21 #######################################################
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
source("model.R")
source("performance_metrics.R")

# Loading ETF returns
returns <- read_csv("data_directory/etfs_rtn.csv")[-1]

# Creating auxiliary matrices and list
We <- 252 # Window size for GARCH estimation
Wt <- nrow(returns) # Total size of window
K <- 1000 # Number of simulations 

# Construct portfolio and benchmarks
mixture_portfolio <- RollingWindowEstimation(returns = returns,
                                             We = 252,
                                             Wt = Wt,
                                             K = K)

gaussian_portfolio <- RollingWindowEstimation(returns = returns,
                                              We = 252,
                                              Wt = Wt,
                                              K = K)

naive_portfolio <- NaiveDiversification(returns)


# Convert the portfolio_returns matrix to an xts object
mixture_portfolio_xts <- xts::xts(mixture_portfolio[,-1], 
                                  order.by = mixture_portfolio$date)

gaussian_portfolio_xts <- xts::xts(gaussian_portfolio[,-1], 
                                   order.by = gaussian_portfolio$date)

naive_portfolio_xts <- xts::xts(naive_portfolio[(We+1):Wt, -1],
                                order.by = naive_portfolio[(We+1):Wt,]$date)