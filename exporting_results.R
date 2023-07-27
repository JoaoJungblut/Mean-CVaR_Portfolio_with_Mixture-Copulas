# Function to save a LaTeX table of statistics summary
SaveSummaryStats <- function(df, filename = "tables/summary_stats_table.txt") {
  
  # SaveSummaryStats: Function to save a LaTeX table of statistics summary.
  # Inputs:
  #   df: A data frame containing asset returns with a 'date' column and individual asset columns.
  #   filename: File path to save the LaTeX table (default = "tables/summary_stats_table.txt").
  
  # Summary statistics for each ETF's log returns
  df[is.na(df)] <- 0
  summary_stats <- data.frame("Symbol" = names(df[-1]),
                              "Min." = colMins(df[-1]),
                              "1st Quantile" = colQuantiles(df[-1], prob = 0.25),
                              "Median." = colQuantiles(df[-1], prob = 0.5),
                              "Mean" = colMeans(df[-1]),
                              "3rd Quantile" = colQuantiles(df[-1], prob = 0.75),
                              "Max." = colMaxs(df[-1]),
                              "Std. Dev" = colSds(df[-1]),
                              "Skewness" = colSkewness(df[-1]),
                              "Kurtosis" = colKurtosis(df[-1]))[-1]
  
  # Convert the data frame to a LaTeX table using xtable
  summary_table <- xtable(summary_stats,
                          caption = "Summary Statistics of ETF Returns",
                          align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c"),
                          digits = c(0, 2, 4, 4, 4, 4, 2, 4, 4, 4))
  
  # Save the LaTeX table as .txt file
  write(as.character(summary_table), file = filename)
  
  cat("Performance graph saved to:", filename, "\n")
}



SavePerformanceTable <- function(all_performance, filename = "tables/performance_table.txt") {
  # SavePerformanceTable: Function to save a LaTeX table of portfolio performance metrics.
  # Inputs:
  #   all_performance: A list containing performance metrics for different portfolio strategies.
  #   filename: File path to save the LaTeX table (default = "tables/performance_table.txt").
  
  # Convert the data frame to a LaTeX table using xtable
  performance_table <- xtable(all_performance,
                              caption = "Portfolio Performance Metrics",
                              align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"),
                              digits = c(2, 4, 4, 4, 4, 4, 4, 4, 4))
  
  # Save the LaTeX table as .txt file
  write(as.character(performance_table), file = filename)
  
  cat("Performance table saved to:", filename, "\n")
}



SaveGraphReturns <- function(df, filename = "figures/returns_figure.png") {
  
  # SaveGraphReturns: Function to save a combined plot of ETF returns in a single figure.
  # Inputs:
  #   df: A data frame containing ETF returns with a 'date' column and individual asset columns.
  #   filename: File path to save the plot as an image (default = "figures/returns_figure.png").
  
  # Convert the data from wide to long format using 'gather' 
  df[is.na(df)] <- 0
  df_long <- tidyr::gather(df, key = "Symbol", value = "Returns", -date)
  
  colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
              "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
              "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d")
  names(colors) <- unique(df_long$Symbol)
  
  # Create the plot using 'ggplot2' and 'facet_wrap'
  p <- ggplot(df_long, aes(x = date, y = Returns, color = Symbol)) +
    geom_line() +
    labs(title = "ETF Returns",
         x = "Date",
         y = "Returns",
         color = "Symbol") +
    theme_bw() +  # White background
    facet_wrap(~ Symbol, scales = "free_y", ncol = 4) +  # Adjust the 'ncol' value as needed
    scale_color_manual(values = colors) +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 10))  # Decrease the axis numbers font size
  
  # Save the plot as a file
  ggsave(filename, plot = p)  # Adjust 'width' and 'height' as needed
  
  cat("Performance graph saved to:", filename, "\n")
}



SavePerformanceGraph <- function(data, filename = "figures/performance_graph.png") {
  
  # SavePerformanceGraph: Function to save a performance summary chart as an image.
  # Inputs:
  #   data: A data frame containing performance metrics for the portfolio.
  #   filename: File path to save the performance graph as an image (default = "figures/performance_graph.png").
  
  # Create the PerformanceSummary chart
  chart <- PerformanceAnalytics::charts.PerformanceSummary(data, 
                                                           main = "Performance chart")
  
  # Save the chart to a file
  dev.copy(png, filename)
  dev.off()
  
  cat("Performance graph saved to:", filename, "\n")
}
