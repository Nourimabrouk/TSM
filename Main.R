"
Time Series Models
Assignment 1

Authors:
Zeus Paraguas
Bart
Jari
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
sv <- read.delim(here('Data', 'sv.dat'))

data <- sv %>% as_tibble() %>% rename(x = X...Pound.Dollar.daily.exchange.rates..sections.9.6.and.14.4)
data

