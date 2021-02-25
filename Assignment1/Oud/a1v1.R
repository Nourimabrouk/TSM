rm(list=ls())

# Imports ----------
library(tidyverse)
library(here)
library(ggplot2)
library(ggthemes)
options(warn=-1)

# Data import--------

setwd(here())
data = Nile

# Functions --------
# Smoothed state

kalman_filter <- function(data, theta, sig_eps, sig_eta){
  y <- as.matrix(data)
  n <- nrow(y)
  a <- rep(0,n)   # filter estimator
  P <- rep(0,n)   # variance of filtered estimator
  v <- rep(0,n)   # prediction error
  F <- rep(0,n)   # variance of prediction error
  K <- rep(0,n)   # kalman gain
  a_y <- rep(0,n) # one-step ahead estimator
  P_y <- rep(0,n) # variance of one-step ahead estimator
  a[1] <- theta[1]
  P[1] <- theta[2]
  
  for (i in 1:n) {
    v[i] <- y[i] - a[i]
    F[i] <- P[i] + sig_eps
    
    if(is.nan(y[i]) || is.na(y[i])){
      K[i] <- 0
      a_y[i] <- a[i]
      P_y[i] <- P[i]
      a[i+1] <- a[i]
      P[i+1] <- P[i] + sig_eta  
      
      print(a[i])
      
    } else{
      
      K[i] <- P[i]/F[i]
      a_y[i] <- a[i] + K[i]*v[i]
      P_y[i]  = P[i]*(1-K[i])
      
      if(i < (n-1)){
        a[i+1] <- a[i] + K[i]*v[i]
        P[i+1] <- P[i] * (1-K[i]) + sig_eta  
      }
    }
  }
  
  # Upper and lower bounds of a
  a_lb <- a-1.645*sqrt(P)
  a_ub <- a+1.645*sqrt(P)
  
  kalman <- data.frame(a, P, v, F, K, a_y, P_y, a_lb, a_ub)
  return(kalman) 
}
smoothed_state <- function(a, P, v, F, K){
  n <- length(v)
  alpha <- rep(0,n) # smoothed stated
  N <- rep(0,n)     # smoothed state error variance
  r <- rep(0,n)
  V <- rep(0,n)     # smoothed state variance
  L <- 1 - K
  
  N[n] <- 0
  r[n] <- 0
  
  for (j in (n-1):2){ #reversed loop
    if (j>1){
      N[j-1] <- (1/F[j]) + L[j]^2 * N[j]
    }
    
    if (is.nan(v[j]) || is.na(v[j])){
      r[j-1] <- r[j]
      
      V[j] <- P[j] - P[j]^2 * N[j-1]   
      alpha[j] <- a[j]+P[j]*r[j-1]
    }
    else {                   
      r[j-1] <- (v[j]/F[j])+L[j]*r[j]
      
      V[j] <- P[j] - P[j]^2 * N[j-1]   
      alpha[j] <- a[j]+P[j]*r[j-1]
    }  
  }
  
  #upper and lower bounds of alpha
  alpha_lb <- alpha - 1.645*sqrt(V)     
  alpha_ub <- alpha + 1.645*sqrt(V)
  
  SmoothedState_df <- data.frame(alpha, N, r, V, alpha_lb, alpha_ub)
  return (SmoothedState_df)
}
disturbances_smoothing <- function(F, v, K, r, N){
  u <- 1/F*V-K*r
  eps <- sig_eps*u
  D <- 1/F + K^2*N
  var_eps <- sig_eps - sig_eps^2*D
  eta <- sig_eta*r
  var_eta <- sig_eta - sig_eta^2*N
  sd_eps <- sqrt(var_eps)
  sd_eta <- sqrt(var_eta)
  disturbance <- data.frame(eps,eta,sd_eps,sd_eta)
  return(disturbance)
}

# -------- Main --------

# Kalman Filter
sig_eps <- 15099
sig_eta <- 1469.1
a_ini <- 0 
P_ini <- 10^7
theta <- c(a_ini, P_ini)

df_kalman_filtered_state <- kalman_filter(data, theta, sig_eps, sig_eta)

# Smoothed State
a <- df_kalman_filtered_state$a
P <- df_kalman_filtered_state$P
v <- df_kalman_filtered_state$v
F <- df_kalman_filtered_state$F
K <- df_kalman_filtered_state$K
df_smoothed_state <- smoothed_state(a, P, v, F, K)

# Disturbances Smoothing
F <- df_kalman_filtered_state$F
V <- df_kalman_filtered_state$v
k <- df_kalman_filtered_state$K
r <- df_smoothed_state$r
N <- df_smoothed_state$N
df_disturbance <- disturbances_smoothing(F, V, K, r, N)

# Results
summary(df_disturbance)
summary(df_smoothed_state)
summary(df_kalman_filtered_state)

# Missing data
missing_values_index <- c(21:40, 61:80)
data_missing <- data
data_missing[missing_values_index] <- NA
data_missing

df_kalman_missing_data <- kalman_filter(data_missing, theta, sig_eps, sig_eta)

a_missing <- df_kalman_missing_data$a
P_missing <- df_kalman_missing_data$P
F_missing <- df_kalman_missing_data$F
v_missing <- df_kalman_missing_data$v
k_missing <- df_kalman_missing_data$K
r_missing <- df_kalman_missing_data$r
N_missing <- df_kalman_missing_data$N

df_smoothed_state_missing_data <- smoothed_state(a_missing, P_missing, v_missing, F_missing, k_missing)

F_missing <- df_kalman_missing_data$F
V_missing <- df_kalman_missing_data$v
k_missing <- df_kalman_missing_data$K
r_missing <- df_smoothed_state_missing_data$r
N_missing <- df_smoothed_state_missing_data$N

df_disturbance_missing_data <- disturbances_smoothing(F_missing, V_missing, k_missing, r_missing, N_missing)


# --------Plotting-------
# Plot KF (2.1)
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))

filtered_state <- df_kalman_filtered_state$a[2:99]
filtered_state_lower_bound <- df_kalman_filtered_state$a_lb[2:99]
filtered_state_upper_bound <- df_kalman_filtered_state$a_ub[2:99]

filtered_variance <- df_kalman_filtered_state$P[2:99]
state_error <- df_kalman_filtered_state$v[2:99]
state_error_variance <- df_kalman_filtered_state$F[2:99]

y_lim_one <- c(min(data), max(data))
y_lim_two <- c(min(filtered_variance), max(filtered_variance))
y_lim_three <- c(min(state_error), max(state_error))
y_lim_four <- c(min(state_error_variance), max(state_error_variance))

plot(ts(filtered_state,start=c(1871, 1)), plot.type="single", ylab="", main="i", ylim=y_lim_one)
lines(ts(filtered_state_lower_bound, start=c(1871, 1)), col="red")
lines(ts(filtered_state_upper_bound, start=c(1871, 1)), col="red")
points(ts(data,start=c(1871, 1)))

plot(ts(filtered_variance, start=c(1871, 1)), plot.type="single", ylab="", main="ii", ylim=y_lim_two)
plot(ts(state_error, start=c(1871, 1)), plot.type="single", ylab="", main="iii", ylim=y_lim_three)
abline(h=0,col="red")
plot(ts(state_error_variance, start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=y_lim_four)

# Plot SS (2.2)

smooth_state <- df_smoothed_state$alpha[2:99]
smooth_state_lower_bound <- df_smoothed_state$alpha_lb[2:99]
smooth_state_upper_bound <- df_smoothed_state$alpha_ub[2:99]

smooth_variance <- df_smoothed_state$V[2:99]
state_error <- df_smoothed_state$r[1:100]
state_error_variance <- df_smoothed_state$N[1:100]

y_lim_one <- c(min(data), max(data))
y_lim_two <- c(min(smooth_variance), max(smooth_variance))
y_lim_three <- c(min(state_error), max(state_error))
y_lim_four <- c(min(state_error_variance), max(state_error_variance))


par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(ts(smooth_state, start=c(1871, 1)), plot.type="single", ylab="", main="i", ylim=y_lim_one)
lines(ts(smooth_state_lower_bound, start=c(1871, 1)), col="red")
lines(ts(smooth_state_upper_bound, start=c(1871, 1)), col="red")
points(ts(data,start=c(1871, 1)), col="red")
plot(ts(df_smoothed_state$V[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="ii", ylim=y_lim_two)
plot(ts(df_smoothed_state$r[1:100], start=c(1871, 1)), plot.type="single", ylab="", main="iii", ylim=y_lim_three)
abline(h=0,col="red")
plot(ts(df_smoothed_state$N[1:100], start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=y_lim_four)

# Plot DS (2.3) 
### Figures 2.3 (ii),(iv) plot standard deviations instead of variances
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))





plot(ts(df_disturbance$eps[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="i", ylim=c(-300,300))
abline(h=0,col="red")
plot(ts(df_disturbance$sd_eps[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="ii", ylim=c(45,65))
plot(ts(df_disturbance$eta[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="iii", ylim=c(-40,40))
abline(h=0,col="red")
plot(ts(df_disturbance$sd_eta[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=c(35,40))
# SKIP PLOT 2.4 - NOT NECESSARY

# Plot MV 2.5
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(ts(df_kalman_missing_data$a[2:99], start=c(1871, 1)), plot.type ="single", ylab="", main="i", ylim=c(500,1400))
plot(ts(df_kalman_missing_data$F[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=c(min(df_kalman_missing_data$F[2:99]),max(df_kalman_missing_data$F[2:99])))
plot(ts(df_smoothed_state_missing_data$alpha[2:99], start=c(1871, 1)), plot.type ="single", ylab="", main="i", ylim=c(500,1400))
plot(ts(df_smoothed_state_missing_data$V[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=c(min(df_kalman_missing_data$F[2:99]),max(df_kalman_missing_data$F[2:99])))

# Plot Forecasting 2.6




# Plot Diagnostic Plots prediction errors 2.7




# Plot Diagnostic Plots auxilliary residuals 2.8



