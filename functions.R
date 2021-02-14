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
create_ylim <- function(vector){
  return(c(min(vector), max(vector)))
}
makeTS <- function(vector){
  ts <- ts(vector, start=c(1871,1))
  return(ts)
}

plotOne <- function(df){
  # Input: df_kalman_filtered_state
  n <- length(data)
  
  filtered_state <- df$a[2:n]
  filtered_state_lb <- df$a_lb[2:n]
  filtered_state_ub <- df$a_ub[2:n]
  
  filtered_variance <- df$P[2:n]
  state_error <- df$v[2:n]
  state_error_variance <- df$F[2:n]
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(filtered_state), plot.type="single", ylab="", main="i", ylim=create_ylim(data))
  lines(makeTS(filtered_state_lb), col="red")
  lines(makeTS(filtered_state_ub), col="red")
  points(makeTS(data))
  
  plot(makeTS(filtered_variance), plot.type="single", ylab="", main="ii", ylim=create_ylim(filtered_variance))
  plot(makeTS(state_error), plot.type="single", ylab="", main="iii", ylim=create_ylim(state_error))
  abline(h=0,col="red")
  plot(makeTS(state_error_variance), plot.type="single", ylab="", main="iv", ylim=create_ylim(state_error_variance))
}
plotTwo <- function(df){
  # input: df_smoothed_state
  # Plot SS (2.2)
  
  # ii) en iv) lijken in de vroegste observaties niet overeen te komen met het boek (pg 22 / 45)
  smooth_state <- df$alpha[2:99]
  smooth_state_lb <- df$alpha_lb[2:99]
  smooth_state_ub <- df$alpha_ub[2:99]
  
  smooth_variance <- df$V[2:99]
  state_error <- df$r[1:100]
  state_error_variance <- df$N[1:100]
  ts(smooth_state, start=c(1871, 1))
  makeTS(smooth_state)
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(makeTS(smooth_state), plot.type="single", ylab="", main="i", ylim=create_ylim(data))
  lines(makeTS(smooth_state_lb), col="red")
  lines(makeTS(smooth_state_ub), col="red")
  points(makeTS(data), col="red")
  plot(makeTS(df$V[2:99]), plot.type="single", ylab="", main="ii", ylim=create_ylim(smooth_variance))
  plot(makeTS(df$r[1:100]), plot.type="single", ylab="", main="iii", ylim=create_ylim(state_error))
  abline(h=0,col="red")
  plot(makeTS(df$N[1:100]), plot.type="single", ylab="", main="iv", ylim=create_ylim(state_error_variance))
  
}
plotThree <- function(df){
  # Input: df_disturbance
  ### Figures 2.3 (ii),(iv) plot standard deviations instead of variances
  ### Missing first observation? Compare with book (pg 25/48)
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(df$eps[2:99]), plot.type="single", ylab="", main="i", ylim=c(-300,300))
  abline(h=0,col="red")
  plot(makeTS(df$sd_eps[2:99]), plot.type="single", ylab="", main="ii", ylim=c(45,65))
  plot(makeTS(df$eta[2:99]), plot.type="single", ylab="", main="iii", ylim=c(-40,40))
  abline(h=0,col="red")
  plot(makeTS(df$sd_eta[2:99]), plot.type="single", ylab="", main="iv", ylim=c(35,40))
}
plotFive <- function(df_k, df_s){
  #Input: df_kalman_missing_data
  # Plot MV 2.5
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(makeTS(df_k$a[2:99]), plot.type ="single", ylab="", main="i", ylim=c(500,1400))
  plot(makeTS(df_k$F[2:99]), plot.type="single", ylab="", main="iv", ylim=create_ylim(df_kalman_missing_data$F[2:99]))
  plot(makeTS(df_s$alpha[2:99]), plot.type ="single", ylab="", main="i", ylim=c(500,1400))
  plot(makeTS(df_s$V[2:99]), plot.type="single", ylab="", main="iv", ylim=create_ylim(df_kalman_missing_data$F[2:99]))

  }
plotSix <- function(df){
  # Plot Forecasting 2.6
  
}
plotSeven <- function(df){
  # Plot Diagnostic Plots prediction errors 2.7
  
}
plotEight <- function(df){
  # Plot Diagnostic Plots auxilliary residuals 2.8
  
}
