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
setwd(here())
# Imports ----------
library(here)
source("functions.R")
source("plotting.R")
library(tidyverse)

options(warn=-1)

# Data import--------

df_flights = read_csv(here('Data', 'total-number-of-flights.csv')) %>% 
  select(2) %>% slice(100:199)
ts_Nile = Nile
# 2.1 Kalman Filter
# Initialises values
# Applies Kalman Filter
# Plots results (figure 2.1)

sig_eps <- 15099
sig_eta <- 1469.1
a_ini <- 0 
P_ini <- 10^7
theta <- c(a_ini, P_ini)

df_kalman_filtered_state <- kalman_filter(ts_Nile, theta, sig_eps, sig_eta)
plotOne(df_kalman_filtered_state)

# 2.2  Smoothed State 
# Creates figure 2.2
df_smoothed_state <- smoothed_state(ts_Nile, df_kalman_filtered_state)
plotTwo(df_smoothed_state)

# 2.3 Disturbance Smoothing
# Creates figure 2.3

df_disturbance <- disturbances_smoothing(df_kalman_filtered_state, df_smoothed_state)
plotThree(df_disturbance)

# 2.5 Missing data
# Creates missing data
# Runs Kalman filter, state smoother, and disturbance smoother
# Creates figure 2.5
missing_values_index <- c(21:40, 61:80)
df_data_missing <- ts_Nile
df_data_missing[missing_values_index] <- NA

df_kalman_missing_data <- kalman_filter(df_data_missing, theta, sig_eps, sig_eta)
df_smoothed_state_missing_data <- smoothed_state(df_data_missing, df_kalman_missing_data)
df_disturbance_missing_data <- disturbances_smoothing(df_kalman_missing_data, df_smoothed_state_missing_data)

plotFive(df_data_missing, df_kalman_missing_data, df_smoothed_state_missing_data)

# 2.6
# Forecasting
# Creates figure 2.6
n_steps <- 30
df_forecasts <- one_step_forecasting(df_kalman_filtered_state, n_steps)
plotSix(ts_Nile, df_kalman_filtered_state, df_forecasts)

# 2.7

# Standardised prediction errors
# Creates figure 2.7

df_predictionerrors <- prediction_errors(df_kalman_filtered_state$v, df_kalman_filtered_state$F)
plotSeven(df_predictionerrors)


# 2.8
# Standardized smoothed residuals
# Creates figure 2.8
df_st_residuals <- stand_smooth_residuals(df_kalman_filtered_state$F,df_kalman_filtered_state$v,
                                          df_kalman_filtered_state$K, df_smoothed_state$r,
                                          df_smoothed_state$N)

plotEight(df_st_residuals,df_predictionerrors)  

# Param optimization
phi_ini <- 0
parameter_estimation <- kalman_parameter_optimizer(df_data, phi_ini)
q_hat <- exp(parameter_estimation$par)

kalman_star <- optimal_kalman_filter(df_data, q_hat)
sig_eps_hat <- calcualte_sig_eps_hat(kalman_star)
sig_eta_hat <- q_hat*sig_eps_hat

theta_hat <- c(q_hat, sig_eps_hat, sig_eta_hat)

cat("The parameter estimates are:")
round(theta_hat, 4)

cat("The log-likelihood value is:")
parameter_estimation$value

cat("iterations:")
parameter_estimation$counts[1]

cat("Exit flag:")
parameter_estimation$convergence # zero indicates succesfull optimization
