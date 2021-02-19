
plotOne <- function(df_data, df_kf){
  "
  Goal: Plot figure 2.1
  Input: df_kalman_filtered_state
  Output: Plot 2.1
  
  "
  # Input: 
  png("1.png", width = 1000, height = 500)
  
  n <- nrow(df_kf)
  
  y <- as.matrix(df_data)[2:n]
  filtered_state <- df_kf$a[2:n]
  filtered_state_lb <- df_kf$a_lb[2:n]
  filtered_state_ub <- df_kf$a_ub[2:n]
  
  filtered_variance <- df_kf$P[2:n]
  state_error <- df_kf$v[2:n]
  state_error_variance <- df_kf$F[2:n]
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(filtered_state, 2), lot.type="single", ylab="", xlab="", main = "(i)",font.main=1, cex.main=1,  ylim=create_ylim(y))
  lines(makeTS(filtered_state_lb, 2), col="red")
  lines(makeTS(filtered_state_ub, 2), col="red")
  points(makeTS(y, 2), pch=20)
  
  plot(makeTS(filtered_variance, 2), plot.type="single", ylab="", xlab="", main = "(ii)",font.main=1, cex.main=1,  ylim=create_ylim(filtered_variance))
  plot(makeTS(state_error, 2), plot.type="single", ylab="", xlab="Days", main="(iii)",font.main=1, cex.main=1,  ylim=create_ylim(state_error))
  abline(h=0,col="red")
  plot(makeTS(state_error_variance, 2), plot.type="single", ylab="", xlab="Days", main = "(iv)",font.main=1, cex.main=1,  ylim=create_ylim(state_error_variance))
  
  dev.off()
  
}

plotTwo <- function(df_data, df_sm){
  "
  Goal: Plot figure 2.4
  Input: df_smoothed_state
  Output: Plot 2.4
  "
  png("2.png", width = 1000, height = 500)
  
  n <- nrow(df_sm)
  
  y <- as.matrix(df_data)[2:n]
  smooth_state <- df_sm$alpha[2:n]
  smooth_state_lb <- df_sm$alpha_lb[2:n]
  smooth_state_ub <- df_sm$alpha_ub[2:n]
  
  smooth_variance <- df_sm$V[2:n]
  state_error <- df_sm$r[2:(n-1)]
  state_error_variance <- df_sm$N[2:(n-1)]
  
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(makeTS(smooth_state, 2), plot.type="single", ylab="", xlab="", main = "(i)",font.main=1, cex.main=1,  ylim=create_ylim(y))
  lines(makeTS(smooth_state_lb, 2), col="red")
  lines(makeTS(smooth_state_ub, 2), col="red")
  points(makeTS(y, 2), pch=20)
  
  plot(makeTS(smooth_variance, 2), plot.type="single", ylab="", xlab="", main = "(ii)",font.main=1, cex.main=1,  ylim=create_ylim(smooth_variance[2:n]))
  plot(makeTS(state_error, 2), plot.type="single", ylab="", xlab="Days", main = "(iii)",font.main=1, cex.main=1,  ylim=create_ylim(state_error))
  abline(h=0,col="red")
  plot(makeTS(state_error_variance, 2), plot.type="single", ylab="", xlab="Days", main = "(iv)",font.main=1, cex.main=1,  ylim=create_ylim(state_error_variance))  
  dev.off()
  
}

plotThree <- function(df){
  "
  Goal: Plot figure 2.3
  Input: df_disturbance
  Output: Plot 2.3 
  "
  png("3.png", width = 1000, height = 500)
  
  n <- nrow(df)
  observation_error <- df$eps[1:n]
  observation_error_variance <- df$sd_eps[1:n]
  state_error <- df$eta[1:n]
  state_error_variance <- df$sd_eta[1:n]     
  
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(makeTS(observation_error, 2), plot.type="single", ylab="", xlab="", main = "(i)",font.main=1, cex.main=1,  ylim=create_ylim(observation_error))
  abline(h=0, col="red")
  plot(makeTS(observation_error_variance, 2), plot.type="single", ylab="", xlab="", main = "(ii)",font.main=1, cex.main=1,  ylim=create_ylim(observation_error_variance))
  plot(makeTS(state_error, 2), plot.type="single", ylab="", xlab="Days", main = "(iii)",font.main=1, cex.main=1,  ylim=create_ylim(state_error))
  abline(h=0, col="red")
  plot(makeTS(state_error_variance, 2), plot.type="single", ylab="", xlab="Days", main = "(iv)",font.main=1, cex.main=1,  ylim=create_ylim(state_error_variance))
  dev.off()
}

plotFive <- function(df_data, df_k, df_s){
  "
  Goal: Plot figure 2.5
  Input: df_data, df_kalman_missing_data, df_smoothed_missing_data
  Output: Plot 2.5
  "
  png("5.png", width = 1000, height = 500)
  
  n <- length(df_data)
  
  filtered_state <- df_k$a[2:n]
  filtered_variance <- df_k$P[2:n]
  
  smoothed_state <- df_s$alpha[2:n]
  smoothed_state_variance <- df_s$V[2:n]
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(makeTS(filtered_state, 2), col="red", plot.type="single", ylab="", xlab="", main = "(i)",font.main=1, cex.main=1,  ylim=create_ylim(df_data))
  lines(makeTS(df_data, 2))
  
  plot(makeTS(filtered_variance, 2), plot.type="single", ylab="", xlab="", main = "(ii)",font.main=1, cex.main=1,  ylim=create_ylim(filtered_variance))
  
  plot(makeTS(smoothed_state, 2), col="red", plot.type="single", ylab="", xlab="Days", main = "(iii)", font.main=1, cex.main=1,  ylim=create_ylim(df_data))
  lines(makeTS(df_data, 2))
  
  plot(makeTS(smoothed_state_variance, 2), plot.type="single", ylab="", xlab="Days", main = "(iv)", font.main=1, cex.main=1,  ylim=create_ylim(smoothed_state_variance))
  dev.off()
  
}

plotSix <- function(df_data, df_filtered, df_forecasts){
  png("6.png", width = 1000, height = 500)
  
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
  
  plot(makeTS(forecast_state, 2), plot.type="single", ylab="", xlab="", main = "(i)",font.main=1, cex.main=1,  ylim=c(30000,200000))
  
  lines(ts(forecast_state_lb, start=c(101,1)), col="red")
  lines(ts(forecast_state_ub, start=c(101,1)), col="red")
  points(makeTS(df_data, 2), pch=20)
  
  plot(makeTS(forecast_variance, 2), plot.type="single", ylab="", xlab="",main = "(ii)",font.main=1, cex.main=1,  ylim=create_ylim(forecast_variance))
  plot(makeTS(forecast_observation, 2), plot.type="single", ylab="", xlab="Days", main = "(iii)",font.main=1, cex.main=1,  ylim=create_ylim(forecast_observation))
  plot(makeTS(forecast_error_variance, 2), plot.type="single", ylab="", xlab="Days", main = "(iv)",font.main=1, cex.main=1,  ylim=create_ylim(forecast_error_variance)) 
  dev.off()
  
}

plotSeven <- function(df){
  "
  Goal: Plot Diagnostic Plots prediction errors 2.7
  Input: 
  Output: Plot 2.7
  "
  png("7.png", width = 1000, height = 500)
  
  n <- nrow(df)
  
  stv <- df[,1]
  temp <- makeTS(stv, 2)
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  
  plot(temp,ylab="", xlab= "Days", main = "(i)",font.main=1, cex.main=1)
  abline(h=0)
  hist(stv, prob=T, col = "white", main = "(ii)",font.main=1, cex.main=1,  xlab ="",ylab="", ylim=c(0,0.45))
  lines(density(stv[-1]))
  qqnorm(c(stv),main = "(iii)",font.main=1, cex.main=1,  pch = 1)
  qqline(c(stv))
  
  plot(acf(stv, lag.max = 10, plot=F), xlim=c(1,10), ylim=c(-1, 2),ci = 0)
  title("(iv)", cex.main = .75, font.main =1)
  dev.off()
  
}
plotEight <- function(df_st_residuals, df_7output){
  "
  Goal: Plot Diagnostic Plots auxilliary residuals 2.8
  Input: 
  Output: Plot 2.8
  "
  png("8.png", width = 1000, height = 500)
  length <- nrow(df_st_residuals)
  
  stv <- df_predictionerrors[-1,1]
  u_star <- makeTS(df_st_residuals[-1,1], 2)
  r_star <- makeTS(df_st_residuals[-1,2], 2)
  
  par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
  plot(u_star,xlab="",ylab="",
       main = "(i)",font.main=1, cex.main=1)
  abline(h=0)
  
  hist(df_st_residuals[-1,1], prob=T, bins = 25, col = "white" ,main = "(ii)",font.main=1, cex.main=1, xlab="",ylab="", na.rm=TRUE)
  lines(density(u_star, na.rm=TRUE), na.rm=TRUE)
  
  plot(r_star,ylab="",xlab="Days" ,main = "(iii)",font.main=1, cex.main=1)
  abline(h=0)
  hist(df_st_residuals[,2], col = "white", prob=T, main = "(iv)",font.main=1, cex.main=1, 
       ylim=c(0,1.2),xlab="Theoretical Quantiles",ylab="", na.rm=TRUE)
  lines(density(r_star, na.rm=TRUE), na.rm=TRUE)
  dev.off()
  
}

