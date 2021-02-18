kalman_filter <- function(data, theta, sig_eps, sig_eta){
  "
  Goal: Apply the kalman filter recursion 
  Input: data, theta, sig_eps, sig_eta
  Output: DF 100x9 - data.frame(a, P, v, F, K, a_y, P_y, a_lb, a_ub)
  
  "
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
      
      
      if(i < (n-1)){
        a[i+1] <- a[i]
        P[i+1] <- P[i] + sig_eta 
      }
      
      a[i+1] <- a[i]
      P[i+1] <- P[i] + sig_eta  
      
      
    } else{
      
      K[i] <- P[i]/F[i]
      a_y[i] <- a[i] + K[i]*v[i]
      P_y[i]  = P[i]*(1-K[i])
      
      if(i < (n-1)){
        a[i+1] <- a[i] + K[i]*v[i]
        P[i+1] <- P[i] * (1-K[i]) + sig_eta  
      }
    }
    
    a[n] <- a[n-1] + K[n-1]*v[n-1]
    P[n] <- P[n-1] * (1-K[n-1]) + sig_eta  
    
    
  }
  
  # Upper and lower bounds of a
  a_lb <- a-1.645*sqrt(P)
  a_ub <- a+1.645*sqrt(P)
  
  kalman <- data.frame(a, P, v, F, K, a_y, P_y, a_lb, a_ub)
  return(kalman) 
}

smoothed_state <- function(df){
  "
  Goal: Compute smoothed state through reverse loop
  Input: df_kalman_filtered_state (output of kalman filter function)
  Output: DF 100x6 - data.frame(alpha, N, r, V, alpha_lb, alpha_ub)
  
  "
  a <- df$a
  P <- df$P
  v <- df$v
  F <- df$F
  K <- df$K
  
  
  n <- length(v)
  alpha <- rep(0,n) # smoothed stated
  N <- rep(0,n)     # smoothed state error variance
  r <- rep(0,n)
  V <- rep(0,n)     # smoothed state variance
  L <- 1 - K
  
  V[1] = P[1]
  
  N[n] <- 0
  r[n] <- 0
  
  
  for (j in n:2){ #reversed loop
    if (j > 1){
      N[j-1] <- (1/F[j]) + L[j]^2 * N[j]
      V[j] <- P[j] - P[j]^2 * N[j-1] 
    }
    if (is.nan(v[j]) || is.na(v[j])){
      r[j-1] <- r[j]
      alpha[j] <- a[j] + P[j]*r[j-1]
    }
    else {                   
      r[j-1] <- (v[j]/F[j])+L[j]*r[j]
      alpha[j] <- a[j]+P[j]*r[j-1]
    }  
  }
  
  #upper and lower bounds of alpha
  alpha_lb <- alpha - 1.645*sqrt(V)     
  alpha_ub <- alpha + 1.645*sqrt(V)
  
  SmoothedState_df <- data.frame(alpha, N, r, V, alpha_lb, alpha_ub)
  return (SmoothedState_df)
}

disturbances_smoothing <- function(dfKalman, dfSmoothed){
  
  "
  Goal: Apply disturbance smoothing
  Input: df_kalman_filtered_state,df_smoothed_state (output of kalman, smoothedstate functions)
  Output: DF 100x4 - data.frame(eps,eta,sd_eps,sd_eta)
  
  
  "
  F <- dfKalman$F
  V <- dfKalman$v
  K <- dfKalman$K
  r <- dfSmoothed$r
  N <- dfSmoothed$N
  
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

create_ylim <- function(vector){
  "
  Goal: Create min-max range from data
  Input: Vector
  Output: List c(min, max)
  "
  return(c(min(vector, na.rm=TRUE), max(vector, na.rm=TRUE)))
}

makeTS <- function(vector,c){
  "
  Goal: Convert vector into time series starting in 1871
  Input: Vector
  Output: Time series
  "
  if(c==1){
    ts <- ts(vector, start=c(1871,1))
  } else {
    ts <-ts(vector,start=c(1970,1))
  }
  return(ts)
}

one_step_ahead_forecast <- function(data, theta, sig_eps, sig_eta, n_steps){
  return(data)
}

forecasting <- function(dfkalman){
  steps <- 30
  n <- nrow(dfkalman)
  yearWithForecast <- seq(n+steps)
  
  a_forecast <- rep(0,steps)
  P_forecast <- rep(0,steps)
  F_forecast <- rep(0,steps)
  
  a_forecast[1] <- array(dfkalman$a)[n]
  P_forecast[1] <- array(dfkalman$P)[n] + sig_eta
  F_forecast[1] <- P_forecast[1] + sig_eps
  
  for (j in 1:(steps-1) ) { 
    a_forecast[j+1] <- a_forecast[j]
    P_forecast[j+1] <- P_forecast[j] + sig_eta
    F_forecast[j+1] <- P_forecast[j+1] + sig_eps
  }
  a_lb_forecast <- a_forecast - 0.675*sqrt(F_forecast)
  a_ub_forecast <- a_forecast + 0.675*sqrt(F_forecast)
  
  forecast <- data.frame(a_forecast, P_forecast, F_forecast, a_lb_forecast,a_ub_forecast)
  return(forecast)
}

prediction_errors <- function(v,F){
  n <- length(v)
  sfe <- rep(0,n) # Standardised forecast error
  
  for (i in 1: n) {
    sfe[i] <- v[i]/sqrt(F[i])
  }
  st_error <- data.frame(sfe,v,F)
  return(st_error)
}

stand_smooth_residuals <- function(F,v,K,r,N){
  u <- 1/F*v-K*r
  D <- 1/F+K^2*N
  u_star <-u/sqrt(D)
  r_star <- r/sqrt(N)
  
  st_residuals <- data.frame(u_star,r_star)
  return(st_residuals)
}

plotOne <- function(df){
  "
  Goal: Plot figure 2.1
  Input: df_kalman_filtered_state
  Output: Plot 2.1
  
  "
  # Input: 
  n <- nrow(df)
  
  filtered_state <- df$a[2:n]
  filtered_state_lb <- df$a_lb[2:n]
  filtered_state_ub <- df$a_ub[2:n]
  
  filtered_variance <- df$P[2:n]
  state_error <- df$v[2:n]
  state_error_variance <- df$F[2:n]
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(filtered_state,1), plot.type="single", ylab="", main="i", ylim=create_ylim(data))
  lines(makeTS(filtered_state_lb,1), col="red")
  lines(makeTS(filtered_state_ub,1), col="red")
  points(makeTS(data,1))
  
  plot(makeTS(filtered_variance,1), plot.type="single", ylab="", main="ii", ylim=create_ylim(filtered_variance))
  plot(makeTS(state_error,1), plot.type="single", ylab="", main="iii", ylim=create_ylim(state_error))
  abline(h=0,col="red")
  plot(makeTS(state_error_variance,1), plot.type="single", ylab="", main="iv", ylim=create_ylim(state_error_variance))
}

plotTwo <- function(df){
  "
  Goal: Plot figure 2.4
  Input: df_smoothed_state
  Output: Plot 2.4
  "
  
  n <- nrow(df)
  smooth_state <- df$alpha[2:n]
  smooth_state_lb <- df$alpha_lb[2:n]
  smooth_state_ub <- df$alpha_ub[2:n]
  
  smooth_variance <- df$V[2:n]
  state_error <- df$r[1:(n-1)]
  state_error_variance <- df$N[1:(n-1)]
  
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(makeTS(smooth_state,1), plot.type="single", ylab="", main="i", ylim=create_ylim(data))
  lines(makeTS(smooth_state_lb,1), col="red")
  lines(makeTS(smooth_state_ub,1), col="red")
  points(makeTS(data,1), col="red")
  plot(makeTS(smooth_variance,1), plot.type="single", ylab="", main="ii", ylim=create_ylim(smooth_variance))
  plot(makeTS(state_error,1), plot.type="single", ylab="", main="iii", ylim=create_ylim(state_error))
  abline(h=0,col="red")
  plot(makeTS(state_error_variance,1), plot.type="single", ylab="", main="iv", ylim=create_ylim(state_error_variance))  
}

plotThree <- function(df){
  "
  Goal: Plot figure 2.3
  Input: df_disturbance
  Output: Plot 2.3 
  "
  n <- nrow(df)
  observation_error <- df$eps[2:n]
  observation_error_variance <- df$sd_eps[2:n]
  state_error <- df$eta[2:n]
  state_error_variance <- df$sd_eta[2:n]     
  
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(makeTS(observation_error,1), plot.type="single", ylab="", main="i", ylim=create_ylim(observation_error))
  abline(h=0, col="red")
  plot(makeTS(observation_error_variance,1), plot.type="single", ylab="", main="ii", ylim=create_ylim(observation_error_variance))
  plot(makeTS(state_error,1), plot.type="single", ylab="", main="iii", ylim=create_ylim(state_error))
  abline(h=0, col="red")
  plot(makeTS(state_error_variance,1), plot.type="single", ylab="", main="iv", ylim=create_ylim(state_error_variance))}
plotFive <- function(df_data, df_k, df_s){
  "
  Goal: Plot figure 2.5
  Input: df_data, df_kalman_missing_data, df_smoothed_missing_data
  Output: Plot 2.5
  "
  n <- length(df_data)
  
  filtered_state <- df_k$a[2:n]
  filtered_variance <- df_k$P[2:n]
  
  smoothed_state <- df_s$alpha[2:n]
  smoothed_state_variance <- df_s$V[2:n]
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(filtered_state,1), col="red", plot.type="single", ylab="", main="i", ylim=create_ylim(df_data))
  lines(makeTS(df_data,1))
  
  plot(makeTS(filtered_variance,1), plot.type="single", ylab="", main="ii", ylim=create_ylim(filtered_variance))
  
  plot(makeTS(smoothed_state,1), col="red", plot.type="single", ylab="", main="iii", ylim=create_ylim(df_data))
  lines(makeTS(df_data,1))
  
  plot(makeTS(smoothed_state_variance,1), plot.type="single", ylab="", main="iv", ylim=create_ylim(smoothed_state_variance))
}

plotSix <- function(df,dv){
  "
  Goal: Plot Forecasting 2.6
  Input: df_kalman_filtered_state, df_forecasting
  Output: Plot 2.6
  "
  m <- nrow(df)
  n <- nrow(dv)

  forecast_state <- c(df$a[2:m],dv$a_forecast[1:n])
  forecast_state_lb <- dv$a_lb_forecast[1:n]
  forecast_state_ub <- dv$a_ub_forecast[1:n]
  
  forecast_variance <- c(df$P[2:m],dv$P_forecast[1:n])
  forecast_observation <- c(df$a[2:m],dv$a_forecast[1:n])
  forecast_error_variance <- c(df$F[2:m],dv$F_forecast[1:n])
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(forecast_state,1), plot.type="single", ylab="", main="i", ylim=create_ylim(data))
  lines(makeTS(forecast_state_lb,2), col="red")
  lines(makeTS(forecast_state_ub,2), col="red")
  points(makeTS(data,1))
  
  plot(makeTS(forecast_variance,1), plot.type="single", ylab="", main="ii", ylim=create_ylim(forecast_variance))
  plot(makeTS(forecast_observation,1), plot.type="single", ylab="", main="iii", ylim=create_ylim(forecast_observation))
  plot(makeTS(forecast_error_variance,1), plot.type="single", ylab="", main="iv", ylim=create_ylim(forecast_error_variance)) 
}

plotSeven <- function(df){
  "
  Goal: Plot Diagnostic Plots prediction errors 2.7
  Input: 
  Output: Plot 2.7
  "
  
  
}
plotEight <- function(df){
  "
  Goal: Plot Diagnostic Plots auxilliary residuals 2.8
  Input: 
  Output: Plot 2.8
  "
  
  
}

