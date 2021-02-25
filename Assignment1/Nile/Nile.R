"
Time Series Models
Assignment 1

Authors:
Zeus Paraguas
Bart
Jari
Nouri Mabrouk 2623401


Main file to run the analysis
"
rm(list=ls())
# Imports ----------
library(here)
source("functionsNile.R")
source("plottingNile.R")
library(tidyverse)

setwd(here())
options(warn=-1)

# Data import--------

df_flights = read_csv(here('Data', 'total-number-of-flights.csv')) %>% 
  select(2) %>% slice(100:199)
ts_Nile = Nile

# 2.1

sig_eps <- 15099
sig_eta <- 1469.1
a_ini <- 0 
P_ini <- 10^7
theta <- c(a_ini, P_ini)

df_kalman_filtered_state <- kalman_filter(ts_Nile, theta, sig_eps, sig_eta)
plotOne(ts_Nile,df_kalman_filtered_state)

# 2.2
df_smoothed_state <- smoothed_state(ts_Nile, df_kalman_filtered_state)
plotTwo(ts_Nile,df_smoothed_state)

# 2.3
df_disturbance <- disturbances_smoothing(df_kalman_filtered_state, df_smoothed_state)
plotThree(df_disturbance)

# 2.5 
missing_values_index <- c(21:40, 61:80)
df_data_missing <- ts_Nile
df_data_missing[missing_values_index] <- NA

df_kalman_missing_data <- kalman_filter(df_data_missing, theta, sig_eps, sig_eta)
df_smoothed_state_missing_data <- smoothed_state(df_data_missing, df_kalman_missing_data)
df_disturbance_missing_data <- disturbances_smoothing(df_kalman_missing_data, df_smoothed_state_missing_data)

plotFive(df_data_missing, df_kalman_missing_data, df_smoothed_state_missing_data)

# 2.6

n_steps <- 30
df_forecasts <- one_step_forecasting(df_kalman_filtered_state, n_steps)
plotSix(ts_Nile, df_kalman_filtered_state, df_forecasts)

# 2.7


df_predictionerrors <- prediction_errors(df_kalman_filtered_state$v, df_kalman_filtered_state$F)
plotSeven(df_predictionerrors)

# 2.8

df_st_residuals <- stand_smooth_residuals(df_kalman_filtered_state$F,df_kalman_filtered_state$v,
                                          df_kalman_filtered_state$K, df_smoothed_state$r,
                                          df_smoothed_state$N)

plotEight(df_st_residuals, df_predictionerrors)  

