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
library(tidyverse)
library(here)
source("functions.R")

options(warn=-1)

# Data import--------

setwd(here("Data"))
data = Nile

# 2.1 Kalman Filter
# Initialises values
# Applies Kalman Filter
# Plots results (figure 2.1)

sig_eps <- 15099
sig_eta <- 1469.1
a_ini <- 0 
P_ini <- 10^7
theta <- c(a_ini, P_ini)

df_kalman_filtered_state <- kalman_filter(data, theta, sig_eps, sig_eta)
plotOne(df_kalman_filtered_state)


# 2.2  Smoothed State 
# Creates figure 2.2

df_smoothed_state <- smoothed_state(df_kalman_filtered_state)
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
data_missing <- data
data_missing[missing_values_index] <- NA

df_kalman_missing_data <- kalman_filter(data_missing, theta, sig_eps, sig_eta)
df_smoothed_state_missing_data <- smoothed_state(df_kalman_missing_data)
df_disturbance_missing_data <- disturbances_smoothing(df_kalman_missing_data, df_smoothed_state_missing_data)

plotFive(data_missing, df_kalman_missing_data, df_smoothed_state_missing_data)
# 2.6
# 2.7
# 2.8
# 2.9
  
