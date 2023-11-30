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
source("garch_estimate.R")
source("copula_estimate.R")
source("portfolio_optimization.R")
source("pipeline_module.R")


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


# Calculate cumulative returns
results <- Pipeline(Ret_inSample, Ret_outofSample, Update, 
                    copulas = c("Gumbel", "t"), K = 1000,
                    Alpha = 0.05, TargetReturn = 0, NumAssets = 8)
cumulative_returns <- cumsum(results)


# Plot the cumulative returns
plot(cumulative_returns$V1, type = "l", col = "black", lwd = 1,
     main = "Portfolio Performance",
     xlab = "Time Period", ylab = "Cumulative Returns")

