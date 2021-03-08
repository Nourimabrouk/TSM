"
Time Series Models
Assignment 2

Authors:
Zeus Paraguas 2624650
Bart
Jari Verbeek 2580924
Nouri Mabrouk 2623401


This code applies CH 2 of the book including the figures presented to our own time series
"

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

# Data import--------

data <- read.delim(here('Data', 'sv.dat'))
stonks <- read_csv(here('Data', 'oxfordmanrealizedvolatilityindices.csv'))

returns <- data %>% 
  mutate(index = 1:nrow(data)) %>% 
  relocate(index) %>% 
  as_tsibble(index = index) %>% 
  rename(x = X...Pound.Dollar.daily.exchange.rates..sections.9.6.and.14.4) %>% 
  mutate(demeaned = (x - mean(x))/100,
         transformed = log(demeaned^2)) 

stonkdata <- stonks %>%   
  filter(Symbol == ".SPX" & year(X1) > 2015) %>% 
  select(X1,close_price, rk_parzen) %>% # replace rk_parzen with realized volatility measure of choice
  rename(Date = X1, Close = close_price, RV = rk_parzen) %>%
  mutate(RV = log(RV)
         ) %>% 
  as_tsibble()

y <- diff(log(stonkdata$Close))
x <- log((y - mean(y))^2)

rv <- stonkdata$RV[-1]
stonks_data1 <- cbind(x, rv)
  
returns
stonkdata

source("functions.R")
sig_eps <- (pi^2)/2
mean_u <- -1.27

par_ini <- c(0.1082, 0.40, -0.207, 0.0)
Z <- 1
H <- sig_eps
T <- par_ini[2]
R <- par_ini[1]
Q <- 1
Beta <- 0

c <- par_ini[3]
d <- mean_u

state_space_parameters <- data.frame(Z, H, T, R, Q, Beta, c, d)
ret_trans <- returns$transformed
res <- state_space_parameter_optimizer(ret_trans, par_ini, state_space_parameters)



h <- rep(0,N)
y_simulate <- rep(0,N)


N <- 10000
phi <- 0.99
sigma <- 0.95
omega <- 0.2

epsilon <- rnorm(N)
eta <- rnorm(N,0,sqrt(sigma))

h[1] <- omega

for (t in 1:N){
  y_simulate[t] <- h[t] + epsilon[t]
  h[t+1] <- omega + phi*h[t] + sigma*eta[t]
}



source("functions.R")
res2 <- state_space_parameter_optimizer(y_simulate, par_ini, state_space_parameters)




source("functions.R")
res2 <- state_space_parameter_optimizer(stonks_data1, par_ini, state_space_parameters)





# e)
# Overview of dataset
stonks %>% head
colnames(stonks)
unique(stonks$Symbol)
range(stonks$X1)

# Data prep

#Quickplots
autoplot(returns, demeaned)
autoplot(returns, transformed)

autoplot(stonkdata, Close)
autoplot(stonkdata, RV)
range(returns$demeaned)

# 14.5 (i)
ggplot(returns, aes(index, demeaned))+
  theme_minimal()+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_y_continuous(breaks = seq(-.05,.05,.025))+
  scale_x_continuous(breaks = seq(0,900,100))
# # 14.5 (ii)
# ggplot(returns, aes(index))+
#   theme_minimal()+
#   geom_point()+
#   # + geom_line(smoothed_estimate k+theta)
#   geom_hline(yintercept = 0)+
#   scale_x_continuous(breaks = seq(0,900,100))

# 14.5 (iii)
# ggplot(returns, aes(index, SE_volmeasure))+
#   theme_minimal()+
#   geom_line()
