rm(list=ls())
# Imports ----------

source("functions.R")
source("Plotting.R")

library(here)
library(tidyverse)
library(tsibble)
library(lubridate)
library(scales)
library(ggplot2)
library(ggthemes)
library(fable)
library(tseries)
library(moments)

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

# Descriptive stats
descriptive_a <- descriptive_stats(returns$x)
descriptive_a
descriptive_b <- descriptive_stats(returns$transformed)
descriptive_b

descriptive_e1 <- descriptive_stats(input_stocks[,1])
descriptive_e1 
descriptive_e2 <- descriptive_stats(input_stocks[,-1])
descriptive_e2 

par_ini <- c(0.1082, 0.98, -0.2, 0.9)
initial_parameters <- initialise_parameters_QML(par_ini)

state_space_parameters <- initial_parameters[[1]]
par_ini <- initial_parameters[[2]]

QML_params_returns <- optimize_parameters(input_returns, par_ini[-4], state_space_parameters, TRUE) # (Print_output = TRUE)
outputKalman_returns <- compute_kalmanfilter(input_returns, QML_params_returns, state_space_parameters)
outputSmooth_returns <- compute_smoothed_state(input_returns, QML_params_returns, outputKalman_returns)

xi_sv <- QML_params_returns[3]/(1 - QML_params_returns[2])
h_t <- outputKalman_returns$h_t
plot_returns_input <- returns %>% mutate(
  alpha = outputSmooth_returns$alpha,
  H_filtered = h_t - xi_sv,
  H_smoothed = outputSmooth_returns$alpha - xi_sv)

#e SPX data and Log RV model extension
par_ini <- c(0.1082, 0.98, -0.2, 0.5)
QML_params_stocks <- optimize_parameters(input_stocks[,1], par_ini[-4], state_space_parameters, TRUE)
QML_params_stocks_rv <- optimize_parameters(input_stocks, par_ini, state_space_parameters, TRUE)

outputKalman_stocks <- compute_kalmanfilter(input_stocks[,1], QML_params_stocks, state_space_parameters)
outputSmooth_stocks <- compute_smoothed_state(input_stocks[,1], QML_params_stocks, outputKalman_stocks)

outputKalman_stocks_rv <- compute_kalmanfilter(input_stocks, QML_params_stocks_rv, state_space_parameters)
outputSmooth_stocks_rv <- compute_smoothed_state(input_stocks, QML_params_stocks_rv, outputKalman_stocks_rv)

#plot 1 
plot(ts(outputSmooth_returns$alpha) , col="red", plot.type="single", ylab="", main="h_t", ylim=c(min(returns$transformed), max(returns$transformed)))
lines(ts(outputSmooth_stocks_rv$alpha))
points(input_stocks[,1], col="black")

log_RV <- input_stocks[,-1]
Beta_hat <- QML_params_stocks_rv[4]


xi_stocks <- QML_params_stocks[3]/(1 - QML_params_stocks[2])
xi_stocks_rv <- QML_params_stocks_rv[3]/(1 - QML_params_stocks_rv[2])
h_t_stock <- outputKalman_stocks$h_t
h_t_stock_rv <- outputKalman_stocks_rv$h_t
H_filtered_stock <- h_t_stock - xi_stocks
H_filtered_stock_rv <- h_t_stock_rv - xi_stocks_rv
H_smoothed <- outputSmooth_stocks$alpha - xi_stocks 

H_smoothed_rv <- outputSmooth_stocks_rv$alpha - xi_stocks_rv 

plot(ts(H_filtered_stock), col="red", plot.type="single", ylab="", main="H_t Filtered")
plot(ts(H_filtered_stock_rv))

plot(ts(H_smoothed), col="red", plot.type="single", ylab="", main="H_t Smoothed")
plot(ts(H_smoothed_rv))


#f

plot_returns_input <- returns %>% mutate(
  alpha = outputSmooth_returns$alpha,
  H_filtered = h_t - xi_sv,
  H_smoothed = outputSmooth_returns$alpha - xi_sv)


#particle_filtered_stock <- particle_filter(stockdata)
plot_stock_input <- stockdata %>% 
  slice(-1) %>% 
  mutate(
    index = 1:1286,
    x = input_stocks[,1],
    alpha = outputSmooth_stocks$alpha,
    alpha_rv = outputSmooth_stocks_rv$alpha,
    h_t_stock = h_t_stock,
    h_t_stock_rv = h_t_stock_rv,
    H_filtered_stock = H_filtered_stock,
    H_filtered_stock_rv = H_filtered_stock_rv,
    particle_filtered_stock = particle_filtered_stock)

saveRDS(plot_returns_input, file = "plot_returns_input.rds")
saveRDS(plot_stock_input, file = "plot_stock_input.rds")

# e) stockdata
# f) ? 

