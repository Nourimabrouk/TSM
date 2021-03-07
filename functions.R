state_space_parameter_optimizer <- function(df_data, phi_ini){
  
  print(phi_ini)
  results <- optim(par=phi_ini, fn=function(par) - GetloglikGauss(df_data, par), method="BFGS")
  
  theta_hat <- results$par

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

GetloglikGauss<- function(data, theta){
  "
  Goal: Compute Gaussian log likelihood defined in equation 7.2 in DK
  Input: data matrix, parameters of interest vector theta
  Output: loglikihood function value
  
  "
  
  y <- as.matrix(data)
  n <- length(y)
  
  kf_state <- KalmanFilterSV(data, theta)
  
  v <- kf_state$v
  F <- kf_state$F
  
  log_density <- -(1/2)*log(2*pi) - (1/2)*log(abs(F)) - (1/2)*(v^2)/F
  log_density[is.nan(log_density)] <- 0
  loglikelihood <- (sum(log_density))
  
  return(loglikelihood)
  
}

KalmanFilterSV <- function(data, theta){
  "
  Goal: Perform Kalman filter for state space model as defined in eq 4.2 and slide 21 of Week III
  Input: data matrix, parameters of interest vector theta
  Output: DF nx6 - data.frame(alpha, N, r, V, alpha_lb, alpha_ub)
  
  "
  X <- as.matrix(data)
  n <- dim(X)[1]
  m <- dim(X)[2]

  y <- X[,1]

  if (m > 1){
    x <- X[,-1]
    
    sig_eta <- theta[1]
    phi <- theta[2]
    omega <- theta[3]
    beta <- theta[4]
    
  } else{
    x <- rep(0, n)
    
    sig_eta <- theta[1]
    phi <- theta[2]
    omega <- theta[3]
    beta <- 0
  }
  
  sig_u <- pi^2/2    
  mean_u <- -1.27
  
  R_t <- 1
  H_t <- sig_u
  Q_t <- sig_eta
  T_t <- phi
  
  
  # Define Kalman filtering matrices
  h <- rep(0, n)
  P <- rep(0, n)
  v <- rep(0, n)
  F <- rep(0, n)
  K <- rep(0, n)
  
  # Define initial values for unconditional mean and variacne resp.
  h[1] <- omega/(1 - T_t) 
  P[1] <- Q_t/(1 - T_t^2) 
  
  for (t in 1:n) {
    #print(cbind(x[t], h[t], y[t]))
    v[t] <- (y[t] - mean_u) - h[t] - x[t]*beta
    F[t] <- P[t] + H_t
    K[t] <- T_t*(P[t]/F[t])

    if(t < n-1){
      h[t+1] <- omega + T_t*h[t] + K[t]*v[t]
      P[t+1] <- T_t^2*P[t] + (Q_t^2)*R_t - K[t]^2*F[t]
    } 
  }
  
  kalmanfiltersv <- data.frame(h, P, v, F, K)
  
  return(kalmanfiltersv)
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

# Particle filter
# ( Ask karim )
particlefilter <- function(data){
  tao = 1000 # burn in
  n = length(data)
  for (t in (tao + 1): n){
    # Discard start up sample / burnin
    
    # Sample a_t (12.16)
    # Compute weights & Normalise -> as 12.17
    # Compute x^_t
    # Resample
    print("particle filter")
    
  }
  }

