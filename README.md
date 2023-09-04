# Exploring Worst-Case Mixture of Copulas for Optimal Mean-CVaR Portfolio Allocation

This repository contains the code and documentation for my final course project at Universidade Federal do Rio Grande do Sul (UFRGS) in economic science. The project focuses on the usage of mixture copulas for portfolio optimization using the worst-case Conditional Value at Risk (CVaR). The analysis is conducted on a dataset of ETF (Exchange-Traded Fund) returns.

## Project Structure

- **data_directory/**: This directory contains the following CSV files:
  - *etfs_rtn.csv*: Contains the returns data of various ETFs used in the analysis.
- **PortofolioRiskOpt/**: This directory serves as the main RStudio project directory and contains the following files:
  - **main.R**: This script implements the main functionality of the project, including data preprocessing, GARCH (Generalized Autoregressive Conditional Heteroskedasticity) modeling for volatility estimation, copula modeling, portfolio optimization, portfolio analysis, and performance metrics calculation.
  - **data_preprocessing.R**: This script contains functions for preprocessing the raw data, including reading and cleaning the ETF returns data.
  - **garch_estimate.R**: This script contains functions for estimating GARCH models to model the volatility of the ETF returns.
  - **copula_estimate.R**: This script contains functions for estimating mixture copulas for the given dataset, which are used for capturing the joint distribution of the ETF returns.
  - **portfolio_optimization.R**: This script contains functions for implementing the portfolio optimization algorithm using the worst-case CVaR approach, aiming to find an optimal allocation of assets that minimizes CVaR.
  - **portfolio_analysis.R**: This script contains functions for estimating portfolios weights in a rolling window.
  - **performance_metrics.R**: This script contains a function for calculating various performance metrics for the optimized portfolios, helping to evaluate the efficiency and risk of the portfolio strategies.
  - **exporting_results.R**: This script contains functions to export the results of the analysis, including summary statistics, performance tables, and graphs.


## Data Collection

The data for this project was collected from the *Historical Market Data - Stooq*. Unfortunately, the data collection code is not included in this repository, but the resulting dataset is available in the *etfs_rtn.csv* file in the `data_directory` directory.

## How to Run

To run this project, ensure that you have R and RStudio installed on your system. Follow these steps:

1. Clone this repository to your local machine.
2. Open the `PortfolioRiskOpt.Rproj` file in RStudio.
3. Install the required R packages if they are not already installed. You can do this by running `install.packages("package_name")` in the R console.
4. Open the `main.R` file and run the script to execute the entire project workflow.
5. The script will perform data preprocessing, GARCH modeling, copula modeling, portfolio optimization, and calculate performance metrics for various portfolio strategies.
6. The results, including summary statistics, performance tables, and graphs, will be saved in the respective directories as specified in the script.

Please note that the execution order of the modules is important, and it is already defined in the `main.R` script. Make sure to follow this order to get meaningful results and analysis.
