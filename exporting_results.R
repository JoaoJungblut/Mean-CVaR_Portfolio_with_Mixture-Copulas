# Function to save a LaTeX table of statistics summary
SaveSummaryStats <- function(df, filename = "tables/summary_stats_table.txt") {
  
  # SaveSummaryStats: Function to save a LaTeX table of statistics summary.
  # Inputs:
  #   df: A data frame containing asset returns with a 'date' column and individual asset columns.
  #   filename: File path to save the LaTeX table (default = "tables/summary_stats_table.txt").
  
  library(xtable)
  
  # Summary statistics for each ETF's log returns
  summary_stats <- summary(df[-1]) # Exclude the date column from summary
  skewness_values <- apply(df[-1], 2, skewness) # Calculate skewness for each column
  kurtosis_values <- apply(df[-1], 2, kurtosis) # Calculate kurtosis for each column
  
  # Create a data frame to store the summary statistics
  summary_df <- data.frame(Metric = c("Mean", "Std. Dev.", "Min", "1st Qu.", "Median", "3rd Qu.", "Max", "Skewness", "Kurtosis"),
                           t(summary_stats), skewness_values, kurtosis_values)
  
  # Convert the data frame to a LaTeX table using xtable
  summary_table <- xtable(summary_df,
                          caption = "Summary Statistics of ETF Returns",
                          align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"),
                          digits = c(2, 4, 4, 4))
  
  # Save the LaTeX table as .txt file
  write(as.character(summary_table), file = filename)
  
  cat("Performance graph saved to:", filename, "\n")
}



SavePerformanceTable <- function(returns, filename = "tables/performance_table.txt") {
  
  # SavePerformanceTable: Function to save a LaTeX table of portfolio performance metrics.
  # Inputs:
  #   returns: A data frame containing asset returns with a 'date' column and individual asset columns.
  #   filename: File path to save the LaTeX table (default = "tables/performance_table.txt").
  
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
  
  cat("Performance graph saved to:", filename, "\n")
}



SaveGraphReturns <- function(df, filename = "figures/returns_figure.png") {
  
  # SaveGraphReturns: Function to save a combined plot of ETF returns in a single figure.
  # Inputs:
  #   df: A data frame containing ETF returns with a 'date' column and individual asset columns.
  #   filename: File path to save the plot as an image (default = "figures/returns_figure.png").
  
  # Convert the data from wide to long format using 'gather' 
  df_long <- tidyr::gather(df, key = "Symbol", value = "Returns", -date)
  
  colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  names(colors) <- unique(df_long$ETF)
  
  # Create the plot using 'ggplot2' and 'facet_wrap'
  p <- ggplot(df_long, aes(x = date, y = Returns)) +
    geom_line() +
    labs(title = "ETF Returns",
         x = "Date",
         y = "Returns",
         color = "Symbol") +
    theme_minimal() +
    facet_wrap(~ ETF, scales = "free_y", ncol = 4) +  # Adjust the 'ncol' value as needed
    scale_color_manual(values = colors) +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  # Add title to each graph inside the facets
  p <- p + labs(title = NULL) + facet_wrap(~ ETF, scales = "free_y", ncol = 4,
                                           labeller = labeller(ETF = function(x) paste("ETF", x)))
  
  # Save the plot as a file
  ggsave(filename, plot = p, width = 15, height = 10, units = "cm")  # Adjust 'width' and 'height' as needed
  
  cat("Performance graph saved to:", filename, "\n")
}



SavePerformanceGraph <- function(data, filename = "figures/performance_graph.png") {
  
  # SavePerformanceGraph: Function to save a performance summary chart as an image.
  # Inputs:
  #   data: A data frame containing performance metrics for the portfolio.
  #   filename: File path to save the performance graph as an image (default = "figures/performance_graph.png").
  
  # Create the PerformanceSummary chart
  chart <- charts.PerformanceSummary(data)
  
  # Convert the chart to a ggplot object
  gg_chart <- as.ggplot(chart)
  
  # Save the ggplot to a file
  ggsave(filename, gg_chart)
  
  cat("Performance graph saved to:", filename, "\n")
}
