descriptive_stats <- function(data){
  number_observations <- length(data)
  mean <- mean(data)
  var <- var(data)
  min <- min(data)
  max <- max(data)
  skew <- skewness(data)
  kurt <- kurtosis(data)
  
  descriptive <- data.frame(number_observations,mean,var,min,max,skew,kurt)
  return(descriptive)
}
  
get_mle_standard_errors <- function(hessian){

  standard_errors <- sqrt(abs(diag(solve(-hessian))))

  return(standard_errors)
  
}

optimize_parameters <- function(df_data, phi_ini, state_space_matrices, print_output){
  
  results <- optim(par=phi_ini, fn=function(par) - compute_loglikelihood(df_data, par, state_space_matrices), method="BFGS",  hessian = TRUE)
  
  theta_hat <- results$par
  hess <- results$hessian
  standard_errors <- get_mle_standard_errors(hess)

  sigma_star <- exp(theta_hat[1])
  phi_star <- exp(theta_hat[2])/(1 + exp(theta_hat[2]))
  omega_star <- theta_hat[3]
  
  if(length(theta_hat) == 4){
    beta_star <- theta_hat[4]
    theta_star <- c(sigma_star, phi_star, omega_star, beta_star)
  } else{
    theta_star <- c(sigma_star, phi_star, omega_star)
  }
  
  if(print_output == TRUE){
  print_optimizer_output(theta_star, standard_errors, results)
  }
  return(theta_star)
  
}

print_optimizer_output <- function(theta_hat, standard_errors, results){
 
  print("The parameter estimates are:")
  print(round(theta_hat, 4))
  
  print("std. errors")
  print(round(standard_errors, 4))
  
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
  
  if(length(theta) == 4){
    beta_star <- theta[4]
    theta_star <- c(sigma_star, phi_star, omega_star, beta_star)
  } else{
    theta_star <- c(sigma_star, phi_star, omega_star)
  }
    
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
    Beta <- theta[4]
  } else{
    x <- rep(0, n)
    Beta <- 0
  }

  # Extract state space model parameter matrices
  R <- theta[1]
  H <- state_space_matrices$H
  Q <- state_space_matrices$Q
  T <- theta[2]
  Z <- state_space_matrices$Z

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
  
  # Final index of arrays  
  h[n] <- c + T*h_t[n-1]            
  P[n] <- T^2*P_t[n-1] + Q*(R^2)

  h_t[n] <- h[n] + P[n]*Z*v[n]/F[n]
  P_t[n] <- P[n] - (P[n]^2)*(Z^2)/F[n]
    
  }
  
  output_kalmanfilter <- data.frame(h, h_t, P, v, F, K)
  
  return(output_kalmanfilter)
}

transform_data <- function(stockdata, returns){
  
  y <- diff(log(stockdata$Close))
  x <- log((y - mean(y))^2)
  rv <- stockdata$RV[-1]
  
  input_matrix_stocks <- cbind(x, rv)
  input_returns <- returns$transformed
  
  stock_data <- cbind(x, stockdata$RV[-1])
  ret_trans <- returns$transformed
  
  return(list(stock_data, ret_trans))
}

initialise_parameters_QML <- function(par_ini){
  
  sig_eps <- (pi^2)/2 # Given in assignment
  mean_u <- -1.27 # Given in assignment
  
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
  
  return(list(state_space_parameters, par_ini))
}

compute_smoothed_state <- function(data, theta, kf){

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
  L <- phi - K*Z
  
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
perform_particlefilter_routine <- function(n, sigma_eta, phi, sigma, theta_t, y_t){
  # Implementation from 14.5.3 DK
  # With bootstrap filter
  # And pg284 resampling
  
  #Sigma eta, phi from previous optimizations
  theta_t = initialise_theta(n, sigma_eta, phi) #(theta_0)
  for (t in 1:n) {
    # Initialise values <- current, draw_values
    # Generate vector (step i of recursion)
    rnorm(n, phi*theta, sigma_eta^2)
    values = draw_values(n, sigma_eta, phi) # Vector of length n : theta_0 as random sample from normal unconditional distribution of theta
    normalised_weights = compute_normalised_weights(sigma, theta, y_t)
    c(att, ptt) = compute_att_ptt(theta, normalised_weights)
    a_t = resampling(att)
    alpha_t = resampling(att)
  } 
  return(a_t)
}

initialise_theta <- function(n, sigma_eta, phi){
  var = sigma_eta^2 / (1 - phi^2)
  theta_0 = rnorm(n, 0, var)
  return(theta_0)
}
draw_values <- function(theta, sigma_eta, phi){
  n = length(theta)
  var = sigma_eta^2 / (1 - phi^2)
  theta_0 = rnorm(n, 0, var)
  return(theta_0)
}

compute_normalised_weights<- function(theta_t, y_t, x_i){
  
  sigma = sqrt(exp(x_i)) # answer to lucas discussion board
  weights = exp ( -log (2*pi*sigma^2 ) / 2  - theta_t / 2 - ( exp(- theta_t) * y_t^2)) / (2 * sigma ^ 2) # 322 DK (ii)   
  normalised_weights = weights / sum ( weights )
  
  return(normalised_weights)
}

compute_att_ptt<- function(weights, theta){
  a_hat_t_t = sum ( weights * theta_t_i)
  p_hat_t_t = sum ( weights * theta_t_i ^ 2 - a_hat_t_t ^ 2 )
  
  return(cbind(a_hat_t_t, p_hat_t_t, weights))
}


resampling<- function(df_att_ptt_weights){
  # testing placeholders
  #   a_hat_t_t = 1:100
  #   p_hat_t_t = 1:100
  #   weights = 100:1 / sum(100:1)
  #   df_att_ptt_weights <- cbind(a_hat_t_t, p_hat_t_t, weights)
  
  # Resampling from pg 284
  
  att = df_att_ptt_weights[,1]
  weights = df_att_ptt_weights[,3]
  weights = abs(rnorm(100,1,.5))
  resampled = sample_n(input, 100, replace = TRUE, weight = input$weights)
  
  return(resampled)
}
