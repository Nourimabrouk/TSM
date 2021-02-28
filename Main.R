"
Time Series Models
Assignment 1

Authors:
Zeus Paraguas
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
library(lubridate)
library(scales)
library(ggplot2)

setwd(here())
options(warn=-1)

# Data import--------

df_flights <- read_csv(here('Data', 'total-number-of-flights.csv'))
data <- read.delim(here('Data', 'sv.dat'))

## if data in prices
# k <- nrow(data)
# ret <- vector()
# for(i in 2:k){
#   ret[i-1] <- 100*log(data[i,1]/data[i-1,1])
# }

return <- sv %>% as_tibble() %>% rename(x = X...Pound.Dollar.daily.exchange.rates..sections.9.6.and.14.4)

# Demean to avoid taking logs of zeros
mean <- mean(return$x)

ret <- ret-mean
plot(ts(ret)) # Divide ret by 100 to obtain plot page 320 of DK

# Make SV-model linear by transformation 
ret_trans <- log(ret^2)
plot(ts(ret_trans)) # Again divide by 100 (Moeten het anders plotten volgens de assignment)




