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


source("functions.R")
par_ini <- c(0.1, 0.98, -0.201, 0.9)
perform_QML_routine(returns, stockdata, par_ini)

# QML / cde?
y <- diff(log(stonkdata$Close))
x <- log((y - mean(y))^2)

rv <- stonkdata$RV[-1]
stonks_data1 <- cbind(x, rv)
  
sig_eps <- (pi^2)/2
mean_u <- -1.27


Z <- 1
H <- sig_eps
T <- par_ini[2]
R <- par_ini[1]
Q <- 1
Beta <- par_ini[4]

c <- par_ini[3]
d <- mean_u

state_space_parameters <- data.frame(Z, H, T, R, Q, Beta, c, d)
ret_trans <- returns$transformed
res <- optimize_parameters(ret_trans, par_ini, state_space_parameters, TRUE)
res2 <- optimize_parameters(stonks_data1, par_ini, state_space_parameters)

# e)






