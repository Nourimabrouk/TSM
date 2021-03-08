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

perform_QML_routine(returns, stockdata)

# ------
sig_eps <- (pi^2)/2
mean_u <- -1.27
par_ini <- c(0.1082, 0.991, -0.207, 0.0)

state_space_parameters <- data.frame(Z = 1,
                                     H = sig_eps,
                                     T = par_ini[2],
                                     R = par_ini[1],
                                     Q = 1,
                                     c = par_ini[3],
                                     Beta = par_ini[4],
                                     d = mean_u)

y <- diff(log(stockdata$Close))
x <- log((y - mean(y))^2)
 
stock_data <- cbind(x, stockdata$RV[-1])
ret_trans <- returns$transformed

res <- optimize_parameters(ret_trans, par_ini, state_space_parameters,TRUE)
res2 <- optimize_parameters(stock_data, par_ini, state_space_parameters,TRUE)

outputKalman <- compute_kalmanfilter(x, res, state_space_parameters)
outputSmooth <- compute_smoothed_state(x,res,outputKalman)







