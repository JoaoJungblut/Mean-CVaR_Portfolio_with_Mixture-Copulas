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
library(parallel) # Do parallel computations


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


# Construct portfolio and benchmarks
mixture_portfolio_1y <- RollingWindowEstimation(returns = returns,
                                                We = 252,
                                                Wt = nrow(returns),
                                                K = 1000,
                                                Mixture = TRUE)
mixture_portfolio_2y <- RollingWindowEstimation(returns = returns,
                                                We = 504,
                                                Wt = nrow(returns),
                                                K = 1000,
                                                Mixture = TRUE)
mixture_portfolio_5y <- RollingWindowEstimation(returns = returns,
                                                We = 1260,
                                                Wt = nrow(returns),
                                                K = 1000,
                                                Mixture = TRUE)
gaussian_portfolio_1y <- RollingWindowEstimation(returns = returns,
                                                 We = 252,
                                                 Wt = nrow(returns),
                                                 K = 1000,
                                                 Mixture = FALSE)
gaussian_portfolio_2y <- RollingWindowEstimation(returns = returns,
                                                 We = 504,
                                                 Wt = nrow(returns),
                                                 K = 1000,
                                                 Mixture = FALSE)
gaussian_portfolio_5y <- RollingWindowEstimation(returns = returns,
                                                 We = 1260,
                                                 Wt = nrow(returns),
                                                 K = 1000,
                                                 Mixture = FALSE)
naive_portfolio <- NaiveDiversification(returns)


# Saving results
mixture_portfolio_1y %>% 
  write_csv("results/mixture_portfolio_1y.csv")
mixture_portfolio_2y %>% 
  write_csv("results/mixture_portfolio_2y.csv")
mixture_portfolio_5y %>% 
  write_csv("results/mixture_portfolio_5y.csv")
gaussian_portfolio_1y %>% 
  write_csv("results/gaussian_portfolio_1y.csv")
gaussian_portfolio_2y %>% 
  write_csv("results/gaussian_portfolio_2y.csv")
gaussian_portfolio_5y %>% 
  write_csv("results/gaussian_portfolio_5y.csv")
naive_portfolio %>% 
  write_csv("results/naive_portfolio.csv")


# Convert the portfolio_returns matrix to an xts object
mixture_portfolio_1y_xts <- xts::xts(mixture_portfolio_1y[,-1], 
                                     order.by = mixture_portfolio_1y$date)
mixture_portfolio_2y_xts <- xts::xts(mixture_portfolio_2y[,-1], 
                                     order.by = mixture_portfolio_2y$date)
mixture_portfolio_5y_xts <- xts::xts(mixture_portfolio_5_y[,-1], 
                                     order.by = mixture_portfolio_5y$date)
gaussian_portfolio_1y_xts <- xts::xts(gaussian_portfolio_1y[,-1], 
                                      order.by = gaussian_portfolio_1y$date)
gaussian_portfolio_2y_xts <- xts::xts(gaussian_portfolio_2y[,-1], 
                                      order.by = gaussian_portfolio_2y$date)
gaussian_portfolio_5y_xts <- xts::xts(gaussian_portfolio_5y[,-1], 
                                      order.by = gaussian_portfolio_5y$date)
naive_portfolio_xts <- xts::xts(naive_portfolio[(253):Wncol(returns), -1],
                                order.by = naive_portfolio[(We+1):Wt,]$date)


# Compute performance
mixture_portfolio_1y_performance <- ComputePerformance(mixture_portfolio_1y_xts)
mixture_portfolio_2y_performance <- ComputePerformance(mixture_portfolio_2y_xts)
mixture_portfolio_5y_performance <- ComputePerformance(mixture_portfolio_5y_xts)
gaussian_portfolio_1y_performance <- ComputePerformance(gaussian_portfolio_1y_xts)
gaussian_portfolio_2y_performance <- ComputePerformance(gaussian_portfolio_2y_xts)
gaussian_portfolio_5y_performance <- ComputePerformance(gaussian_portfolio_5y_xts)
naive_portfolio_performance <- ComputePerformance(naive_portfolio_xts)


# Merge the three performance results into one list
all_performance <- list(
  mixture_1y = mixture_portfolio_1y_performance,
  mixture_2y = mixture_portfolio_2y_performance,
  mixture_5y = mixture_portfolio_5y_performance,
  gaussian_1y = gaussian_portfolio_1y_performance,
  gaussian_2y = gaussian_portfolio_2y_performance,
  gaussian_5y = gaussian_portfolio_5y_performance,
  naive = naive_portfolio_performance
)

merged_portfolio_1y <- merge.xts(mixture_portfolio_1y_xts,
                                 gaussian_portfolio_1y_xts,
                                 naive_portfolio_xts, join = "inner")
merged_portfolio_2y <- merge.xts(mixture_portfolio_2y_xts,
                                 gaussian_portfolio_2y_xts,
                                 naive_portfolio_xts, join = "inner")
merged_portfolio_5y <- merge.xts(mixture_portfolio_5y_xts,
                                 gaussian_portfolio_5y_xts,
                                 naive_portfolio_xts, join = "inner")


# Saving results
SaveSummaryStats(df = returns, 
                 filename = "tables/etf_summary_stats_table.txt")
SavePerformanceTable(returns = all_performance, 
                     filename = "tables/etf_performance_table.txt")
SaveGraphReturns(df = returns, 
                 filename = "figures/etf_returns_figure.png")
SavePerformanceGraphs(data = merged_portfolio_1y, 
                      filename = "figures/etf_performance_1y_graph.png")
SavePerformanceGraphs(data = merged_portfolio_2y, 
                      filename = "figures/etf_performance_2y_graph.png")
SavePerformanceGraphs(data = merged_portfolio_5y, 
                      filename = "figures/etf_performance_5y_graph.png")

