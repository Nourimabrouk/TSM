source("currPlotting.R")

plotOne(ts_flights[-1,], df_kalman_filtered_state[-1,])
plotTwo(ts_flights, df_smoothed_state)
plotThree(df_disturbance %>% slice(-c(1,2)))
plotFive(df_data_missing[-1,], df_kalman_missing_data %>% slice(-1), df_smoothed_state_missing_data %>% slice(-1))
plotSix(ts_flights[-1,], df_kalman_filtered_state[-1,], df_forecasts[-1,])
plotSeven(df_predictionerrors %>% slice(-c(1:2)))
plotEight(df_st_residuals %>% slice(-1), df_predictionerrors %>% slice(-1))


