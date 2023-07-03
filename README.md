# The usage of dynamic copulas for portfolio optimization, worst-case CVaR approach with cardinality constraints

This repository contains the code and documentation for my final course project at Universidade Federal do Rio Grande do Sul (UFRGS) in economic science. The project focuses on the usage of dynamic copulas for portfolio optimization using the worst-case Conditional Value at Risk (CVaR) approach with cardinality constraints. The analysis is performed on a dataset of S&P 500 returns since the year 2000.

## Project Structure

- **data_directory/**: This directory contains the following CSV files:
  - *adj_close.csv*: Contains the adjusted close prices of the S&P 500 components.
  - *log_rtn.csv*: Contains the logarithmic returns of the S&P 500 components.
  - *missing_symbols.csv*: Contains the symbols of any missing components in the dataset.
  - *sp500_tickers.csv*: Contains the ticker symbols of all S&P 500 components.
- **PortofolioRiskOpt/**: This directory contains the RStudio project file `PortofolioRiskOpt.Rproj`, which serves as the entry point for the project.
  - **main.R**: Implements the main functionality of the project.
  - **data_preprocessing.R**: Contains functions for preprocessing the raw data.
  - **garch_estimate.R**: Contains functions for estimating GARCH (Generalized Autoregressive Conditional Heteroskedasticity) models for volatility estimation.
  - **copula_estimate.R**: Contains functions for estimating dynamic copulas for the given dataset.
  - **portfolio_optimization.R**: Contains functions for implementing the portfolio optimization algorithm using the worst-case CVaR approach with cardinality constraints.
  - **performance_metrics.R**: Contains functions for calculating various performance metrics for the optimized portfolios.
  - **PortfolioRiskOpt_Documentation.Rmd**: This file contains the R Markdown document, which provides a detailed explanation of the code and functions used in the project.

## Data Collection

The data for this project was collected from the Alpha Vantage API using Python. Unfortunately, the data collection code is not included in this repository, but the resulting dataset is available in the *log_rtn.csv* file in the `data_directory` directory.

## How to Run

To run this project, you need to have R and RStudio installed. Follow these steps:

1. Clone the repository to your local machine.
2. Open the `PortfolioRiskOpt.Rproj` file in RStudio.
3. Install the required R packages if they are not already installed. You can do this by running `install.packages("package_name")` in the R console.
4. Open the `PortfolioRiskOpt_Documentation.Rmd` file and knit it to generate the project report.

Please note that the execution order of the modules may be important, so make sure to follow the logical order outlined in the project report.
