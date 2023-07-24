GetReturns <- function(tickers, start_date) {
  
  # GetReturns: Function to import data and calculate log returns for a list of tickers from a specified start date.
  # Inputs:
  #   tickers: A character vector containing the tickers of assets to fetch data for.
  #   start_date: The starting date to fetch data from (YYYY-MM-DD format).
  # Output:
  #   A data frame containing the date and log returns for each asset in separate columns.
  
  # Fetch data for the specified tickers from the start_date
  data <- tidyquant::tq_get(tickers, from = start_date) %>% 
    # Select relevant columns for date, symbol, and adjusted closing prices
    dplyr::select("date", "symbol", "adjusted") %>%
    # Group data by symbol
    dplyr::group_by(symbol) %>% 
    # Calculate log returns for each asset
    dplyr::mutate(return = log(adjusted) - log(dplyr::lag(adjusted))) %>% 
    # Ungroup the data
    dplyr::ungroup() %>% 
    # Drop the 'adjusted' column
    dplyr::select(-adjusted) %>% 
    # Remove rows with NAs
    stats::na.omit() %>% 
    # Spread the data into columns with each asset's log returns
    tidyr::spread(key = symbol, value = return) %>% 
    # Convert all numeric columns to numeric data type (in case of character columns)
    dplyr::mutate(dplyr::across(where(is.numeric)))
  
  return(data)
}


