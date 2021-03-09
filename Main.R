rm(list=ls())
# Imports ----------

source("functions.R")
#source("Plotting.R")

library(here)
library(tidyverse)
library(tsibble)
library(lubridate)
library(scales)
library(ggplot2)
library(fable)

setwd(here())
options(warn=-1)

# Data import & cleaning --------

data <- read.delim(here('Data', 'sv.dat'))
stocks <- read_csv(here('Data', 'oxfordmanrealizedvolatilityindices.csv'))

returns <- data %>% 
  mutate(index = 1:nrow(data)) %>% 
  relocate(index) %>% 
  as_tsibble(index = index) %>% 
  rename(x = X...Pound.Dollar.daily.exchange.rates..sections.9.6.and.14.4) %>% 
  mutate(demeaned = (x - mean(x))/100,
         transformed = log(demeaned^2)) 

stockdata <- stocks %>%   
  filter(Symbol == ".SPX" & year(X1) > 2015) %>% 
  select(X1,close_price, rk_parzen) %>% # replace rk_parzen with realized volatility measure of choice
  rename(Date = X1, Close = close_price, RV = rk_parzen) %>%
  mutate(RV = log(RV)) %>% 
  as_tsibble()
# ab
transformed_df = transform_data(stockdata, returns)
input_stocks = transformed_df[[1]]
input_returns = transformed_df[[2]]

#c / e -> need to separate these?
initial_parameters = initialise_parameters_QML()

state_space_parameters = initial_parameters[[1]]
par_ini = initial_parameters[[2]]

QML_params_returns <- optimize_parameters(input_returns, par_ini, state_space_parameters, TRUE) # (Print_output = TRUE)
QML_params_stocks <- optimize_parameters(input_stocks, par_ini, state_space_parameters, TRUE)

#d
outputKalman_returns <- compute_kalmanfilter(input_returns, QML_params_returns, state_space_parameters)
outputSmooth_returns <- compute_smoothed_state(input_returns, QML_params_returns, outputKalman_returns)

outputKalman_stocks <- compute_kalmanfilter(input_stocks[,1], QML_params_stocks, state_space_parameters)
outputSmooth_stocks <- compute_smoothed_state(input_stocks[,1],QML_params_stocks, outputKalman_stocks)

#f
# From tutorial 5:
n = 100; 
omega = -0.088; phi = 0.991; sigma_eta = 0.084
#perform_particlefilter_routine(n, sigma_eta, phi, sigma, theta_t, y_t)
