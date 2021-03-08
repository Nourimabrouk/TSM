optimize_parameters <- function(df_data, phi_ini, state_space_matrices, print_output){
  
  results <- optim(par=phi_ini, fn=function(par) - compute_loglikelihood(df_data, par, state_space_matrices), method="BFGS")
  
  theta_hat <- results$par

  sigma_star <- exp(theta_hat[1])
  phi_star <- exp(theta_hat[2])/(1 + exp(theta_hat[2]))
  omega_star <- theta_hat[3]
  beta_star <- theta_hat[4]
  
  theta_star <- c(sigma_star, phi_star, omega_star, beta_star)
  
  if(print_output == TRUE){
  print_optimizer_output(theta_star, results)
  }
  return(theta_hat)
  
}

print_optimizer_output <- function(theta_star, results){
 
   print("The parameter estimates are:")
  print(round(theta_star, 4))
  
  print("The log-likelihood value is:")
  print(results$value)
  
  print("iterations:")
  print(results$counts[1])
  
  cat("Exit flag:")
  print(results$convergence) # zero indicates succesfull optimization
  
  # Reporting estimates:
  # Parameter, logtransformed parameter, SE(logtransformedparameter)
  # As in page 320 of text book
  
}

compute_loglikelihood<- function(data, theta, state_space_matrices){

  y <- as.matrix(data)
  n <- length(y)
  
  sigma_star <- exp(theta[1])
  phi_star <- exp(theta[2])/(1 + exp(theta[2]))
  omega_star <- theta[3]
  beta_star <- theta[4]
  
  theta_star <- c(sigma_star, phi_star, omega_star, beta_star)
  
  kf_state <- compute_kalmanfilter(data, theta_star, state_space_matrices)
  
  v <- kf_state$v
  F <- kf_state$F
  
  log_density <- -(1/2)*log(2*pi) - (1/2)*log(abs(F)) - (1/2)*(v^2)/F
  #log_density[is.nan(log_density)] <- 0
  loglikelihood <- (sum(log_density))
  
  return(loglikelihood)
  
}

compute_kalmanfilter <- function(data, theta, state_space_matrices){

  
  X <- as.matrix(data)
  n <- dim(X)[1]
  m <- dim(X)[2]

  y <- X[,1]

  if (m > 1){
    x <- X[,-1]
    
  } else{
    x <- rep(0, n)
    
  }

  # Extract state space model parameter matrices
  R <- theta[1]
  H <- state_space_matrices$H
  Q <- state_space_matrices$Q
  T <- theta[2]
  Z <- state_space_matrices$Z
  Beta <- theta[4]
  
  c <- theta[3]
  d <- state_space_matrices$d
  
  # Define Kalman filtering matrices
  h <- rep(0, n)
  P <- rep(0, n)
  v <- rep(0, n)
  F <- rep(0, n)
  K <- rep(0, n)
  
  h_t <- rep(0, n)
  P_t <- rep(0, n)

  # Define initial values for unconditional mean and variance resp.
  h[1] <- c/(1 - T) 
  P[1] <- Q/(1 - T^2) 
  
  for (t in 1:n) {
    #print(cbind(x[t], h[t], y[t]))
    v[t] <- (y[t] - d) - Z*h[t] - x[t]*Beta
    F[t] <- Z^2*P[t] + H
    K[t] <- T*(P[t]/F[t])

    h_t[t] <- h[t] + P[t]*Z*v[t]/F[t]
    P_t[t] <- P[t] - (P[t]^2)*(Z^2)/F[t]

    if(t < n-1){
      h[t+1] <- c + T*h_t[t]
      P[t+1] <- T^2*P_t[t] + Q*(R^2)
    } 
  }
  
  output_kalmanfilter <- data.frame(h, P, v, F, K)
  
  return(output_kalmanfilter)
}

perform_QML_routine = function(returns, stockdata){
  
  # Create transformed data matrix
  y = diff(log(stockdata$Close))
  x <- log((y - mean(y))^2)
  rv <- stockdata$RV[-1]
  
  input_matrix_stocks <- cbind(x, rv)
  input_returns <- returns$transformed
  
  # Initialise parameters
  sig_eps <- (pi^2)/2 # Given in assignment
  mean_u <- -1.27 # Given in assignment
  
  par_ini <- c(0.1082, 0.991, -0.207, 0.0)
  
  Beta <- 0
  
  state_space_parameters <- data.frame(
    Q = 1,
    Z = 1,
    H = sig_eps,
    R = par_ini[1],
    T = par_ini[2],
    c = par_ini[3],
    d = mean_u,
    Beta = par_ini[4]
  )
  res <- optimize_parameters(input_returns, par_ini, state_space_parameters, TRUE) # (Print_output = TRUE)
  res2 <- optimize_parameters(input_matrix_stocks, par_ini, state_space_parameters, TRUE)
  
  stock_data <- cbind(x, stockdata$RV[-1])
  ret_trans <- returns$transformed
  
  params_returns <- optimize_parameters(input_returns, par_ini, state_space_parameters,TRUE)
  outputKalman_returns <- compute_kalmanfilter(input_returns, params_returns, state_space_parameters)
  outputSmooth_returns <- compute_smoothed_state(input_returns,params_returns, outputKalman)
  
  params_stocks <- optimize_parameters(input_matrix_stocks, par_ini, state_space_parameters,TRUE)
  outputKalman_stocks <- compute_kalmanfilter(input_matrix_stocks[,1], params_returns, state_space_parameters)
  outputSmooth_stocks <- compute_smoothed_state(input_matrix_stocks[,1],params_returns,outputKalman)
  
}

compute_smoothed_state <- function(data, theta, kf){
  "
  Goal: Compute smoothed state through reverse loop 
  Input: theta, kf (Output of KalmanfilterSV)
  Output: DF nx6 - data.frame(alpha, N, r, V, alpha_lb, alpha_ub)
  
  "
  y <- as.matrix(data)
  h <- kf$h
  P <- kf$P
  v <- kf$v
  F <- kf$F
  K <- kf$K
  
  phi <- theta[2]
  Z <- 1
  
  n <- length(v)
  alpha <- rep(0, n) # smoothed stated
  N <- rep(0, n)     # smoothed state error variance
  r <- rep(0, n)
  V <- rep(0, n)     # smoothed state variance
  L <- phi-K*Z
  
  N[n] <- 0
  r[n] <- 0
  
  for (j in (n):2){ #backward recursion
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
    alpha[j] <- h[j] + P[j]*r[j-1]
  }
  
  N[1] <- (1/F[2]) + (L[2]^2) * N[2]
  N_0 <- (1/F[1]) + (L[1]^2) * N[1]
  
  V[1] <- P[1] - (P[1]^2)*N_0 
  
  r_0 <- (v[1]/F[1]) + L[1]*r[1]
  alpha[1] <- h[1] + P[1]*r_0
  
  Smoothedstate <- data.frame(alpha, r, N, V)
  
  return (Smoothedstate)
}
