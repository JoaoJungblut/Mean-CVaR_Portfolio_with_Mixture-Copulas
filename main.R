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
library(rugarch)       # Univariate GARCH modeling
library(fGarch)        # Multivariate GARCH modeling
library(copula)        # Copula modeling
library(Rsolnp)        # Nonlinear optimization
library(fPortfolio)    # Portfolio optimization
library(PerformanceAnalytics) # Performance metrics
library(xts) # Time series object
library(xtable) # Create LaTex tables
library(ggplot2) # Produce graph 


# Importing modules
source("data_preprocessing.R")
source("garch_estimate.R")
source("copula_estimate.R")
source("portfolio_optimization.R")
source("portfolio_analysis.R")
source("performance_metrics.R")
source("exporting_results.R")


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
                                             K = K,
                                             Mixture = TRUE)
gaussian_portfolio <- RollingWindowEstimation(returns = returns,
                                              We = 252,
                                              Wt = Wt,
                                              K = K,
                                              Mixture = FALSE)
naive_portfolio <- NaiveDiversification(returns)


# Convert the portfolio_returns matrix to an xts object
mixture_portfolio_xts <- xts::xts(mixture_portfolio[,-1], 
                                  order.by = mixture_portfolio$date)
gaussian_portfolio_xts <- xts::xts(gaussian_portfolio[,-1], 
                                   order.by = gaussian_portfolio$date)
naive_portfolio_xts <- xts::xts(naive_portfolio[(We+1):Wt, -1],
                                order.by = naive_portfolio[(We+1):Wt,]$date)


# Compute performance
mixture_portfolio_performance <- ComputePerformance(mixture_portfolio_xts)
gaussian_portfolio_performance <- ComputePerformance(gaussian_portfolio_xts)
naive_portfolio_performance <- ComputePerformance(naive_portfolio_xts)


# Merge the three performance results into one list
merged_performance <- list(
  mixture = mixture_portfolio_performance,
  gaussian = gaussian_portfolio_performance,
  naive = naive_portfolio_performance
)


# Saving results
SaveSummaryStats(df = returns, filename = "tables/etf_summary_stats_table.txt")
SavePerformanceTable(returns = portfolios, filename = "tables/etf_performance_table.txt")
SaveGraphReturns(df = returns, filename = "figures/etf_returns_figure.png")
SavePerformanceGraphs(data = portfolios, filename = "figures/etf_performance_graph.png")


# Generate graph
PerformanceAnalytics::charts.PerformanceSummary(mixture_portfolio_xts)
PerformanceAnalytics::charts.PerformanceSummary(gaussian_portfolio_xts)
PerformanceAnalytics::charts.PerformanceSummary(naive_portfolio_xts)