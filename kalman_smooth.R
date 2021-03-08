# QML / cde?
y <- diff(log(stockdata$Close))
x <- log((y - mean(y))^2)

rv <- stockdata$RV[-1]
stonks_data1 <- cbind(x, rv)

sig_eps <- (pi^2)/2
mean_u <- -1.27

par_ini <- c(0.1082, 0.991, -0.207, 0.0)
Z <- 1
H <- sig_eps
T <- par_ini[2]
R <- par_ini[1]
Q <- 1
Beta <- 0

c <- par_ini[3]
d <- mean_u

state_space_parameters <- data.frame(Z, H, T, R, Q, Beta, c, d)
ret_trans <- returns$transformed

res <- optimize_parameters(ret_trans, par_ini, state_space_parameters,TRUE)
res2 <- optimize_parameters(stonks_data1, par_ini, state_space_parameters,TRUE)

outputKalman <- compute_kalmanfilter(x, res, state_space_parameters)
outputSmooth <- compute_smoothed_state(x,res,outputKalman)
