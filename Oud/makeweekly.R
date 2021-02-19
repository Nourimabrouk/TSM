rm(list=ls())
setwd(here())
# Imports ----------
library(here)
source("functions.R")
source("plotting.R")
library(tidyverse)
library(lubridate)
library(ggplot2)
library(zoo)
options(warn=-1)

# Data import--------

ts_flights = read_csv(here('Data', '3yearsflights.csv')) %>% select(2) %>% slice(1:779)

dates = seq(from = as.Date("2019-01-01"), to = as.Date("2021-02-17"), by = 'day')
df_flights <- tibble(dates, ts_flights) %>% transmute(Date = ymd(dates), Flights = Flights)


df_flights %>% as_tibble() %>% View()
colnames(df_flights)

df_weekly <- df_flights %>%
  group_by(year = year(Date), week = week(Date)) %>% summarise(value = sum(Flights)) %>% tail(100) %>% ungroup() %>% mutate(index = 1:100)

df_weekly

ggplot(df_weekly) +
  geom_line(aes(x=index, y=value))
