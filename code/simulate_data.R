## Get simulated data
library(BIPnet)

set.seed(1)
simulation_results <- Simulate(setting=1)
data_list <- list(simulation_results$X1, 
                  simulation_results$X2, 
                  simulation_results$Y)
indic_var <- c(0, 0, 1) 
method <- "BIP" # method without grouping information to start
group_list <- NULL # default when grouping information not included

## Scale data
# TODO: Breakout scaling into separate step & choose whether or not to scale by the Frobenius norm to account for differing n_features
data_list <- lapply(data_list, scale)

## Get results
bip_0 <- BIP(dataList = data_list, IndicVar = indic_var, Method = method)

## Store simulated data and results for development/ testing
today_date <- Sys.Date()
fname_data <- paste0("data/", today_date, "_simulated_data.rds")
saveRDS(data_list, fname_data)
fname_results <- paste0("data/", today_date, "_simulated_results.rds")
saveRDS(bip_0, fname_results)