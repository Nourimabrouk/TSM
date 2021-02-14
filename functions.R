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
  
  N[n] <- 0
  r[n] <- 0
  
  for (j in n:2){ #reversed loop
    if (j > 1){
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
disturbances_smoothing <- function(dfKalman, dfSmoothed){
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
