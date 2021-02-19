
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
  
  plot(makeTS(filtered_state,1), lot.type="single", ylab="", main = "(i)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(df_Nile))
  lines(makeTS(filtered_state_lb,1), col="grey")
  lines(makeTS(filtered_state_ub,1), col="grey")
  points(makeTS(df_Nile, 1), pch=20)
  
  plot(makeTS(filtered_variance,1), plot.type="single", ylab="", main = "(ii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(filtered_variance))
  plot(makeTS(state_error,1), plot.type="single", ylab="", main = "(iii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(state_error))
  abline(h=0,col="grey")
  plot(makeTS(state_error_variance,1), plot.type="single", ylab="", main = "(iv)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(state_error_variance))
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
  
  smooth_variance <- df$V[1:n]
  state_error <- df$r[1:(n-1)]
  state_error_variance <- df$N[1:(n-1)]
  
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(makeTS(smooth_state,1), plot.type="single", ylab="", main = "(i)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(df_Nile))
  lines(makeTS(smooth_state_lb,1), col="grey")
  lines(makeTS(smooth_state_ub,1), col="grey")
  points(makeTS(df_Nile,1), pch=20)
  
  plot(makeTS(smooth_variance,1), plot.type="single", ylab="", main = "(ii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(smooth_variance[2:n]))
  plot(makeTS(state_error,1), plot.type="single", ylab="", main = "(iii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(state_error))
  abline(h=0,col="grey")
  plot(makeTS(state_error_variance,1), plot.type="single", ylab="", main = "(iv)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(state_error_variance))  
}

plotThree <- function(df){
  "
  Goal: Plot figure 2.3
  Input: df_disturbance
  Output: Plot 2.3 
  "
  n <- nrow(df)
  observation_error <- df$eps[1:n]
  observation_error_variance <- df$sd_eps[1:n]
  state_error <- df$eta[1:n]
  state_error_variance <- df$sd_eta[1:n]     
  
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(makeTS(observation_error,1), plot.type="single", ylab="", main = "(i)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(observation_error))
  abline(h=0, col="grey")
  plot(makeTS(observation_error_variance,1), plot.type="single", ylab="", main = "(ii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(observation_error_variance))
  plot(makeTS(state_error,1), plot.type="single", ylab="", main = "(iii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(state_error))
  abline(h=0, col="grey")
  plot(makeTS(state_error_variance,1), plot.type="single", ylab="", main = "(iv)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(state_error_variance))}

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
  smoothed_state_variance <- df_s$V[1:n]
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(filtered_state,1), col="grey", plot.type="single", ylab="", main = "(i)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(df_Nile))
  lines(makeTS(df_Nile,1))
  
  plot(makeTS(filtered_variance,1), plot.type="single", ylab="", main = "(ii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(filtered_variance))
  
  plot(makeTS(smoothed_state,1), col="grey", plot.type="single", ylab="", main = "(iii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(df_Nile))
  lines(makeTS(df_Nile,1))
  
  plot(makeTS(smoothed_state_variance,1), plot.type="single", ylab="", main = "(iv)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(smoothed_state_variance))
}

plotSix <- function(df_data, df_filtered, df_forecasts){
  "
  Goal: Plot Forecasting 2.6
  Input: df_kalman_filtered_state, df_forecasting
  Output: Plot 2.6
  "
  n <- nrow(df_filtered)
  j <- nrow(df_forecasts)
  
  forecast_state <- c(df_filtered$a[2:n], df_forecasts$a_forecast[1:j])
  forecast_state_lb <- df_forecasts$a_lb_forecast[1:j]
  forecast_state_ub <- df_forecasts$a_ub_forecast[1:j]
  
  forecast_variance <- c(df_filtered$P[2:n], df_forecasts$P_forecast[1:j])
  forecast_observation <- c(df_filtered$a[2:n], df_forecasts$a_forecast[1:j])
  forecast_error_variance <- c(df_filtered$F[2:n], df_forecasts$F_forecast[1:j])
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(forecast_state,1), plot.type="single", ylab="", main = "(i)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(df_Nile))
  lines(makeTS(forecast_state_lb,2), col="grey")
  lines(makeTS(forecast_state_ub,2), col="grey")
  points(makeTS(df_Nile, 1), pch=20)
  
  plot(makeTS(forecast_variance, 1), plot.type="single", ylab="", main = "(ii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(forecast_variance))
  plot(makeTS(forecast_observation, 1), plot.type="single", ylab="", main = "(iii)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(forecast_observation))
  plot(makeTS(forecast_error_variance, 1), plot.type="single", ylab="", main = "(iv)",font.main=1, cex.main=.75, adj = 0, ylim=create_ylim(forecast_error_variance)) 
}

plotSeven <- function(df){
  "
  Goal: Plot Diagnostic Plots prediction errors 2.7
  Input: 
  Output: Plot 2.7
  "
  n <- nrow(df)
  
  stv <- df[,1]
  temp <- makeTS(stv,1)
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(temp,xlab="",ylab="", main = "(i)",font.main=1, cex.main=.75, adj = 0)
  abline(h=0)
  hist(stv, prob=T, col = "white", main = "(ii)",font.main=1, cex.main=.75, adj = 0, xlab ="",ylab="")
  lines(density(stv[-1]))
  qqnorm(c(stv),main = "(iii)",font.main=1, cex.main=.75, adj = 0, pch = 1)
  qqline(c(stv))
  
  plot(acf(stv[-c(1,2)],main = "(iv)",font.main=1, cex.main=.75, adj = 0, lag.max = 10,plot=F), xlim=c(1,10), ylim=c(-1,1),type="h", ci = 0)
}
plotEight <- function(df_st_residuals, df_7output){
  "
  Goal: Plot Diagnostic Plots auxilliary residuals 2.8
  Input: 
  Output: Plot 2.8
  "
  length <- nrow(df_st_residuals)
  
  stv <- df_7output[,1]
  u_star <- makeTS(df_st_residuals[,1],1)
  r_star <- makeTS(df_st_residuals[,2],1)
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(u_star,xlab="",ylab="",
       main = "(i)",font.main=1, cex.main=.75, adj = 0)
  abline(h=0)
  hist(df_st_residuals[,1], prob=T, col = "white" ,ylim=c(0,0.5),
       main = "(ii)",font.main=1, cex.main=.75, adj = 0,xlab="",ylab="", na.rm=TRUE)
  lines(density(u_star, na.rm=TRUE), na.rm=TRUE)
  
  plot(r_star ,xlab="",ylab="" ,main = "(iii)",font.main=1, cex.main=.75, adj = 0)
  abline(h=0)
  hist(df_st_residuals[,2], col = "white", prob=T, main = "(iv)",font.main=1, cex.main=.75, adj = 0,
       ylim=c(0,1.2),xlab="",ylab="", na.rm=TRUE)
  lines(density(r_star, na.rm=TRUE), na.rm=TRUE)
}
