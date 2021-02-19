"
Time Series Models
Assignment 1

Authors:
Zeus Paraguas
Bart
Jari
Nouri Mabrouk 2623401


This code applies CH 2 of the book including the figures presented to our own time series
"

rm(list=ls())
setwd(here())
# Imports ----------
library(here)
source("Functions.R")
source("Plotting.R")
library(tidyverse)

options(warn=-1)

# Data import--------

ts_flights <- read_csv(here('Data', 'total-number-of-flights.csv')) %>% 
  select(2) %>% slice(32:131) %>% ts()

phi_ini <- 0
param_hat <- kalman_parameter_optimizer(ts_flights, phi_ini)

sig_eps <- param_hat[2]
sig_eta <- param_hat[3]

a_ini <- 0 
P_ini <- 10^5
theta <- c(a_ini, P_ini)

#2.1
df_kalman_filtered_state <- kalman_filter(ts_flights, theta, sig_eps, sig_eta)
#2.2
df_smoothed_state <- smoothed_state(ts_flights, df_kalman_filtered_state)
#2.3
df_disturbance <- disturbances_smoothing(df_kalman_filtered_state, df_smoothed_state)
#2.5
missing_values_index <- c(21:40, 61:80)
df_data_missing <- ts_flights
df_data_missing[missing_values_index] <- NA

df_kalman_missing_data <- kalman_filter(df_data_missing, theta, sig_eps, sig_eta)
df_smoothed_state_missing_data <- smoothed_state(df_data_missing, df_kalman_missing_data)
df_disturbance_missing_data <- disturbances_smoothing(df_kalman_missing_data, df_smoothed_state_missing_data)
#2.6
n_steps <- 30
df_forecasts <- one_step_forecasting(df_kalman_filtered_state, n_steps)
#2.7
df_predictionerrors <- prediction_errors(df_kalman_filtered_state$v, df_kalman_filtered_state$F)
#2.8
df_st_residuals <- stand_smooth_residuals(df_kalman_filtered_state$F,df_kalman_filtered_state$v,
                                          df_kalman_filtered_state$K, df_smoothed_state$r,
                                          df_smoothed_state$N)

#Plot results
plotOne(ts_flights, df_kalman_filtered_state)
plotTwo(ts_flights, df_smoothed_state)
plotThree(df_disturbance)
plotFive(df_data_missing, df_kalman_missing_data, df_smoothed_state_missing_data)
plotSix(ts_flights, df_kalman_filtered_state, df_forecasts)
plotSeven(df_predictionerrors)
plotEight(df_st_residuals, df_predictionerrors)  

