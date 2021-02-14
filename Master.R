"
TSM Assignment 1

Figures:
2.1 Done
2.2 Done 
2.3 First observation missing? Something not completely right atleast
2.5 Panel (iv) missing
2.6 - 2.9 Convert to packageless implementation
"

# Imports ----------
library(tidyverse)
library(here)
library(ggplot2)
library(ggthemes)
options(warn=-1)

source("functions.R")
source("plotting.R")

# Data import--------

setwd(here())
data = Nile

# Kalman Filter
sig_eps <- 15099
sig_eta <- 1469.1
a_ini <- 0 
P_ini <- 10^7
theta <- c(a_ini, P_ini)

df_kalman_filtered_state <- kalman_filter(data, theta, sig_eps, sig_eta)
# Smoothed State
df_smoothed_state <- smoothed_state(df_kalman_filtered_state)
# Disturbance Smoothing
df_disturbance <- disturbances_smoothing(df_kalman_filtered_state, df_smoothed_state)

# Missing data
missing_values_index <- c(21:40, 61:80)
data_missing <- data
data_missing[missing_values_index] <- NA

df_kalman_missing_data <- kalman_filter(data_missing, theta, sig_eps, sig_eta)
df_smoothed_state_missing_data <- smoothed_state(df_kalman_missing_data)
df_disturbance_missing_data <- disturbances_smoothing(df_kalman_missing_data, df_smoothed_state_missing_data)
