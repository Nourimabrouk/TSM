# In comments: function name in python code

# Initial values
initial_values<- function(){
  
}
#KalmanFilter
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
    
  }====
  
  kalman_star <- data.frame(v[2:n], F_star[2:n])
  
  return(kalman_star)
  
}

# GetOptKalman
kalman_parameter_optimizer <- function(df_data, phi_ini){
  
  results <- optim(par=phi_ini, fn=function(par) - gauss_loglik_dc(par, df_data), method='BFGS')
  
  q_hat <- exp(results$par)
  kalman_star <- optimal_kalman_filter(df_data, q_hat)
  sig_eps_hat <- calculate_sig_eps_hat(kalman_star)
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

# is getIniTheta 
calculate_sig_eps_hat <- function(kalman_star){
  n <- nrow(kalman_star)
  
  F_star <- kalman_star$F_star
  v <- kalman_star$v
  
  sig_eps <- (1/(n - 1))*sum((v^2)/F_star)
  
  return(sig_eps)
}

#GetloglikGauss
GetloglikGauss<- function(){
  
}
# LogLikGauss
gauss_loglik_dc <- function(theta, data){
  y <- as.matrix(data)
  n <- length(y)
  
  phi <- theta
  q <- exp(phi)
  
  kalman_star <- optimal_kalman_filter(data, q)
  F_star <- kalman_star$F_star
  sig_eps <- calculate_sig_eps_hat(kalman_star)
  
  loglik <- -(n/2)*log(2*pi) - ((n-1)/2) - ((n-1)/2)*(log(sig_eps)) - (1/2)*sum(log(F_star))
  
  return(loglik)
  
}

# KalmanFilterSV
GetloglikGauss<- function(){
  
}
# SmoothedState
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

# LoadData
GetloglikGauss<- function(){
  
}