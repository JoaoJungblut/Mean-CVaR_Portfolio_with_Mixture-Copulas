# Function to save performance metrics as a LaTeX table in a .txt file
SavePerformanceTable <- function(returns, filename = "tables/performance_table.txt") {
  # Compute performance ratios using the original function
  performance_ratios <- ComputePerformance(returns)
  
  # Create a data frame to store the performance metrics and corresponding names
  performance_df <- data.frame(
    Metric = c("Annualized Return",
               "Annualized Std. Dev.",
               "Sharpe Ratio (Annualized)",
               "Sortino Ratio",
               "Omega Sharpe Ratio",
               "VaR (97.5%)",
               "CVaR (97.5%)",
               "Semi-Deviation",
               "Worst Drawdown"),
    Value = performance_ratios
  )
  
  # Convert the data frame to a LaTeX table using xtable
  performance_table <- xtable(performance_df,
                              caption = "Portfolio Performance Metrics",
                              align = c("l", "c"),
                              digits = c(2, 4))
  
  # Save the LaTeX table as .txt file
  write(as.character(performance_table), file = filename)
}

# Example usage:
# Assuming you have a variable 'portfolio_returns' containing the portfolio returns data
SavePerformanceTable(portfolio_returns, "tables/performance_table.txt")


# Function to plot ETF returns and save in a single figure
PlotReturns <- function(df, filename = "etf_returns_figure.png") {
  
  # Convert the data from wide to long format using 'gather' 
  df_long <- tidyr::gather(df, key = "ETF", value = "Returns", -"date")
  
  # Create the plot using 'ggplot2' and 'facet_wrap'
  p <- ggplot(df_long, aes(x = {{date_column}}, y = LogReturns)) +
    geom_line() +
    labs(title = "ETF Returns",
         x = "Date",
         y = "Log Returns") +
    theme_minimal() +
    facet_wrap(~ ETF, scales = "free_y", ncol = 4)  # Adjust the 'ncol' value as needed
  
  # Save the plot as a file
  ggsave(filename, plot = p, width = 15, height = 10, units = "cm")  # Adjust 'width' and 'height' as needed
}

