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

source("Functions.R")
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
## if data in prices
# k <- nrow(data)
# ret <- vector()
# for(i in 2:k){
#   ret[i-1] <- 100*log(data[i,1]/data[i-1,1])
# }

returns <- data %>% 
  mutate(index = 1:nrow(data)) %>% 
  relocate(index) %>% 
  as_tsibble(index = index) %>% 
  rename(x = X...Pound.Dollar.daily.exchange.rates..sections.9.6.and.14.4) %>% 
  mutate(demeaned = (x - mean(x)) / 100 ,
         transformed = log(demeaned ^ 2) ) 

stonkdata <- stonks %>%   
  filter(Symbol == ".AEX" & year(X1) > 2013) %>% 
  select(X1,close_price,rk_parzen) %>% # replace rk_parzen with realized volatility measure of choice
  rename(Date = X1, Close = close_price, RV = rk_parzen) %>%
  mutate(RV = log(RV)) %>% 
  as_tsibble()

returns
stonkdata

source("Functions.R")
par_ini <- c(2, 0.8, 0.9)
res <- state_space_parameter_optimizer(ret_trans, par_ini)

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
