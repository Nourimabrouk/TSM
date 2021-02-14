plot <- function(){
  # Plot KF (2.1)
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))

filtered_state <- df_kalman_filtered_state$a[2:99]
filtered_state_lower_bound <- df_kalman_filtered_state$a_lb[2:99]
filtered_state_upper_bound <- df_kalman_filtered_state$a_ub[2:99]

filtered_variance <- df_kalman_filtered_state$P[2:99]
state_error <- df_kalman_filtered_state$v[2:99]
state_error_variance <- df_kalman_filtered_state$F[2:99]

y_lim_one <- c(min(data), max(data))
y_lim_two <- c(min(filtered_variance), max(filtered_variance))
y_lim_three <- c(min(state_error), max(state_error))
y_lim_four <- c(min(state_error_variance), max(state_error_variance))

plot(ts(filtered_state,start=c(1871, 1)), plot.type="single", ylab="", main="i", ylim=y_lim_one)
lines(ts(filtered_state_lower_bound, start=c(1871, 1)), col="red")
lines(ts(filtered_state_upper_bound, start=c(1871, 1)), col="red")
points(ts(data,start=c(1871, 1)))

plot(ts(filtered_variance, start=c(1871, 1)), plot.type="single", ylab="", main="ii", ylim=y_lim_two)
plot(ts(state_error, start=c(1871, 1)), plot.type="single", ylab="", main="iii", ylim=y_lim_three)
abline(h=0,col="red")
plot(ts(state_error_variance, start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=y_lim_four)

# Plot SS (2.2)

# ii) en iv) lijken in de vroegste observaties niet overeen te komen met het boek (pg 22 / 45)
smooth_state <- df_smoothed_state$alpha[2:99]
smooth_state_lower_bound <- df_smoothed_state$alpha_lb[2:99]
smooth_state_upper_bound <- df_smoothed_state$alpha_ub[2:99]

smooth_variance <- df_smoothed_state$V[2:99]
state_error <- df_smoothed_state$r[1:100]
state_error_variance <- df_smoothed_state$N[1:100]

y_lim_one <- c(min(data), max(data))
y_lim_two <- c(min(smooth_variance), max(smooth_variance))
y_lim_three <- c(min(state_error), max(state_error))
y_lim_four <- c(min(state_error_variance), max(state_error_variance))

par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(ts(smooth_state, start=c(1871, 1)), plot.type="single", ylab="", main="i", ylim=y_lim_one)
lines(ts(smooth_state_lower_bound, start=c(1871, 1)), col="red")
lines(ts(smooth_state_upper_bound, start=c(1871, 1)), col="red")
points(ts(data,start=c(1871, 1)), col="red")
plot(ts(df_smoothed_state$V[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="ii", ylim=y_lim_two)
plot(ts(df_smoothed_state$r[1:100], start=c(1871, 1)), plot.type="single", ylab="", main="iii", ylim=y_lim_three)
abline(h=0,col="red")
plot(ts(df_smoothed_state$N[1:100], start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=y_lim_four)

# Plot DS (2.3) 
### Figures 2.3 (ii),(iv) plot standard deviations instead of variances
### Missing first observation? Compare with book (pg 25/48)
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))

plot(ts(df_disturbance$eps[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="i", ylim=c(-300,300))
abline(h=0,col="red")
plot(ts(df_disturbance$sd_eps[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="ii", ylim=c(45,65))
plot(ts(df_disturbance$eta[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="iii", ylim=c(-40,40))
abline(h=0,col="red")
plot(ts(df_disturbance$sd_eta[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=c(35,40))
# SKIP PLOT 2.4 - NOT NECESSARY

# Plot MV 2.5
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(ts(df_kalman_missing_data$a[2:99], start=c(1871, 1)), plot.type ="single", ylab="", main="i", ylim=c(500,1400))
plot(ts(df_kalman_missing_data$F[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=c(min(df_kalman_missing_data$F[2:99]),max(df_kalman_missing_data$F[2:99])))
plot(ts(df_smoothed_state_missing_data$alpha[2:99], start=c(1871, 1)), plot.type ="single", ylab="", main="i", ylim=c(500,1400))
plot(ts(df_smoothed_state_missing_data$V[2:99], start=c(1871, 1)), plot.type="single", ylab="", main="iv", ylim=c(min(df_kalman_missing_data$F[2:99]),max(df_kalman_missing_data$F[2:99])))
}
