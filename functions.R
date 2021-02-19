kalman_filter <- function(data, theta, sig_eps, sig_eta){
  "
  Goal: Apply the kalman filter recursion 
  Input: data, theta, sig_eps, sig_eta
  Output: DF nx9 - data.frame(a, P, v, F, K, a_y, P_y, a_lb, a_ub)
  
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
    F[i] <- P[i] + sig_eps
    v[i] <- y[i] - a[i]
    
    if(is.nan(y[i]) || is.na(y[i])){
      K[i] <- 0
      a_y[i] <- a[i]
      P_y[i] <- P[i]
      
      if(i < (n-1)){
        a[i+1] <- a[i]
        P[i+1] <- P[i] + sig_eta
      }
      
    } else{
      K[i] <- P[i]/F[i]
      a_y[i] <- a[i] + K[i]*v[i]
      P_y[i] <- P[i]*(1 - K[i])
      
      if(i < (n-1)){
        a[i+1] <- a[i] + K[i]*v[i]
        P[i+1] <- P[i]*(1 - K[i]) + sig_eta
      }

    }
    
    a[n] <- a[n-1] + K[n-1]*v[n-1]
    P[n] <- P[n-1]*(1 - K[n-1]) + sig_eta 
    
    
  }
  
  # Upper and lower bounds of a
  a_lb <- a-1.645*sqrt(P)
  a_ub <- a+1.645*sqrt(P)
  
  kalman <- data.frame(a, P, v, F, K, a_y, P_y, a_lb, a_ub)
  
  return(kalman) 
}

smoothed_state <- function(df_data, df_kf){
  "
  Goal: Compute smoothed state through reverse loop
  Input: df_kf (output of kalman filter function)
  Output: DF nx6 - data.frame(alpha, N, r, V, alpha_lb, alpha_ub)
  
  "
  y <- as.matrix(df_data)
  a <- df_kf$a
  P <- df_kf$P
  v <- df_kf$v
  F <- df_kf$F
  K <- df_kf$K
  
  n <- length(v)
  alpha <- rep(0, n) # smoothed stated
  N <- rep(0, n)     # smoothed state error variance
  r <- rep(0, n)
  V <- rep(0, n)     # smoothed state variance
  L <- 1-K

  N[n] <- 0
  r[n] <- 0
  
  for (j in n:2){ #backward recursion
      N[j-1] <- (1/F[j]) + (L[j]^2) * N[j]
      V[j] <- P[j] - (P[j]^2)*N[j-1] 
      
      if (is.nan(y[j]) || is.na(y[j])){
        N[j-1] <- N[j]
        V[j] <- P[j] - (P[j]^2)*N[j-1] 
        
        r[j-1] <- r[j]
      }
      else {                   
        r[j-1] <- (v[j]/F[j]) + L[j]*r[j]
      }
    alpha[j] <- a[j] + P[j]*r[j-1]
  }
  
  N[1] <- (1/F[2]) + (L[2]^2) * N[2]
  N_0 <- (1/F[1]) + (L[1]^2) * N[1]
  
  V[1] <- P[1] - (P[1]^2)*N_0 
  
  r_0 <- (v[1]/F[1]) + L[1]*r[1]
  alpha[1] <- a[1] + P[1]*r_0
  
  
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
  disturbance <- data.frame(eps, eta, sd_eps, sd_eta)
  
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

makeTS <- function(vector, c){
  "
  Goal: Convert vector into time series starting in 1871
  Input: Vector
  Output: Time series
  "
  if(c==1){
    ts <- ts(vector, start=c(1871, 1))
  } else {
    ts <-ts(vector, start=c(1970, 1))
  }
  return(ts)
}

one_step_forecasting <- function(dfkalman, j_steps){

  n <- nrow(dfkalman)
  yearWithForecast <- seq(n + j_steps)
  
  a_forecast <- rep(0, j_steps)
  P_forecast <- rep(0, j_steps)
  F_forecast <- rep(0, j_steps)
  
  a_forecast[1] <- array(dfkalman$a)[n]
  P_forecast[1] <- array(dfkalman$P)[n] + sig_eta
  
  for (j in 1:(j_steps-1)) { 
    F_forecast[j] <- P_forecast[j] + sig_eps
    a_forecast[j+1] <- a_forecast[j]
    P_forecast[j+1] <- P_forecast[j] + sig_eta

  }
  
  F_forecast[j_steps] <- P_forecast[j_steps] + sig_eps

  a_lb_forecast <- a_forecast - 0.675*sqrt(F_forecast)
  a_ub_forecast <- a_forecast + 0.675*sqrt(F_forecast)
  
  forecast <- data.frame(a_forecast, P_forecast, F_forecast, a_lb_forecast, a_ub_forecast)
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

stand_smooth_residuals <- function(F, v, K, r, N){
  u <- 1/F*v-K*r
  D <- 1/F+K^2*N
  u_star <-u/sqrt(D)
  r_star <- r/sqrt(N)
  
  df_st_residuals <- data.frame(u_star, r_star)
  return(df_st_residuals)
}




########### PARAMETER ESTIMATION ###################

kalman_parameter_optimizer <- function(df_data, phi_ini){

  results <- optim(par=phi_ini, fn=function(par) - gauss_loglik_dc(par, df_data), method='BFGS')
  
  q_hat <- exp(results$par)
  kalman_star <- optimal_kalman_filter(df_data, q_hat)
  sig_eps_hat <- calcualte_sig_eps_hat(kalman_star)
  sig_eta_hat <- q_hat*sig_eps_hat
  
  theta_hat <- c(q_hat, sig_eps_hat, sig_eta_hat)
  
  print("The parameter estimates are:")
  print(round(theta_hat, 4))
  
  print("The log-likelihood value is:")
  print(results$value)
  
  print("iterations:")
  print(results$counts[1])
  
  cat("Exit flag:")
  print(results$convergence) # zero indicates succesfull optimization

  return(theta_hat)
  
}


optimal_kalman_filter <- function(df_data, q){
  y <- as.matrix(df_data)
  n <- length(y)
  
  a <- rep(0, n)     # smoothed state error variance
  v <- rep(0, n)
  F_star <- rep(0, n)     # smoothed state variance
  P_star <- rep(0, n)
  K <- rep(0, n)
  
  
  a[2] <- y[1]
  P_star[2] <- q + 1
  
  for(t in 2:n){
    v[t] <- y[t] - a[t]
    F_star[t] <- P_star[t] + 1
    K[t] <- P_star[t]/F_star[t]
    
    a[t+1] <- a[t] + K[t]*v[t]
    P_star[t+1] <- P_star[t]*(1 - K[t]) + q
    
  }

  kalman_star <- data.frame(v[2:n], F_star[2:n])
  
  return(kalman_star)
  
}

calcualte_sig_eps_hat <- function(kalman_star){
  n <- nrow(kalman_star)
  
  F_star <- kalman_star$F_star
  v <- kalman_star$v
  
  sig_eps <- (1/(n - 1))*sum((v^2)/F_star)
  
  return(sig_eps)
}

gauss_loglik_dc <- function(theta, data){
  y <- as.matrix(data)
  n <- length(y)
  
  phi <- theta
  q <- exp(phi)
  
  kalman_star <- optimal_kalman_filter(data, q)
  F_star <- kalman_star$F_star
  sig_eps <- calcualte_sig_eps_hat(kalman_star)
  
  loglik <- -(n/2)*log(2*pi) - ((n-1)/2) - ((n-1)/2)*(log(sig_eps)) - (1/2)*sum(log(F_star))

  return(loglik)
  
}


