# Importing data and calculating returns
GetReturns <- function(tickers, start_date) {
  data <- tidyquant::tq_get(tickers, from = start_date) %>% 
    dplyr::select("date", "symbol", "adjusted") %>%
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(return = log(adjusted) - log(dplyr::lag(adjusted))) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-adjusted) %>% 
    tidyr::spread(key = symbol, value = return) %>% 
    dplyr::mutate(dplyr::across(where(is.numeric), ~replace(., is.na(.), 0)))
  
  return(data)
}


