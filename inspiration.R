###### Earlier code  ------

rm(list = ls())

library(KFAS)

model <-
  SSModel(Nile ~ SSMtrend(1, Q = list(matrix(NA))) , H = matrix(NA))
model <- fitSSM(inits = c(log(var(Nile)) , log(var(Nile))) ,
                model = model,
                method = 'BFGS')$model
out <-
  KFS(
    model,
    filtering = 'state',
    smoothing = c('state', 'disturbance'),
    simplify = F
  )

par(mfrow = c(2, 2))
conf.bands <- cbind(out$a, as.vector(out$a) +
                      (sqrt(cbind(out$P)) * qnorm(0.95)) %*% (t(c(-1, 1))))
temp <- cbind(Nile, conf.bands[-1,])
cols <- c("grey", "blue", "red", "red")
lwd <- c(1, 2, 1, 1)
lty <- c(1, 1, 2, 2)
plot.ts(
  temp,
  plot.type = "single",
  col = cols,
  lwd = lwd,
  lty = lty,
  xlab = "",
  ylab = "",
  main = "i"
)
legend(
  "topright",
  legend = c("Observation Data", "Filtered state", "Confidence intervals"),
  col = c("grey", "blue", "red"),
  lty = c(1, 1, 2),
  lwd = c(1, 2, 1),
  bty = "n",
  cex = 0.9
)
temp <- ts(c(out$P)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "ii"
)
temp <- ts(c(out$v[-1]), start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "iii"
)
abline(h = 0, col = "grey")
temp <- ts(c(out$F)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "iv"
)

par(mfrow = c(2, 2))
conf.bands <- cbind(out$alphahat, as.vector(out$alphahat) +
                      (sqrt(cbind(out$V)) * qnorm(0.95)) %*% (t(c(-1, 1))))
temp <- cbind(Nile, conf.bands)
cols <- c("grey", "blue", "red", "red")
lwd <- c(1, 2, 1, 1)
lty <- c(1, 1, 2, 2)
plot.ts(
  temp,
  plot.type = "single",
  col = cols,
  lwd = lwd,
  lty = lty,
  xlab = "",
  ylab = "",
  main = "i"
)
legend(
  "topright",
  legend = c("Observation Data", "Smoothed state", "Confidence intervals"),
  col = c("grey", "blue", "red"),
  lty = c(1, 1, 2),
  lwd = c(1, 2, 1),
  bty = "n",
  cex = 0.7
)
temp <- ts(c(out$V), start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "ii"
)
temp <- ts(c(out$r)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "iii"
)
abline(h = 0, col = "grey")
temp <- ts(c(out$N)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "iv"
)

par(mfrow = c(2, 2))
temp <- ts(c(out$epshat), start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "i"
)
abline(h = 0, col = "grey")
temp <- ts(c(out$V_eps), start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "ii"
)

temp <- ts(c(out$etahat), start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "iii"
)
abline(h = 0, col = "grey")
temp <- ts(c(out$V_eta), start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "iv"
)

set.seed(1)
n <- length(Nile)
theta.uncond <- numeric(n)
theta.uncond[1] <- Nile[1]
for (i in 2:n) {
  theta.uncond[i] <-
    theta.uncond[i - 1] + rnorm(1, 0, sqrt(out$model$Q))
}
theta.cond <- simulateSSM(model, type = c("states"))
eps.sim <- simulateSSM(model, type = c("epsilon"))
eta.sim <- simulateSSM(model, type = c("eta"))

par(mfrow = c(2, 2))
ylim = c(400, max(c(
  out$alphahathat, theta.cond, theta.uncond
)) + 50)
plot(
  1:100,
  Nile,
  type = "l",
  col = "grey",
  ylim = ylim,
  main = "i"
)
points(1:100,
       out$alphahat ,
       type = "l",
       col = "green",
       lwd = 2)
points(1:100,
       c(theta.uncond),
       type = "l",
       col = "red",
       lwd = 2)
leg <- c("actual data", "smoothed estimate", "unconditional sample")
legend(
  "topleft",
  leg,
  col = c("grey", "green", "red"),
  lwd = c(1, 2, 2),
  cex = 0.7,
  bty = "n"
)
plot(
  1:100,
  Nile,
  type = "l",
  col = "grey",
  ylim = ylim,
  main = "ii"
)
points(1:100,
       out$alphahat ,
       type = "l",
       col = "green",
       lwd = 2)
points(1:100,
       c(theta.cond),
       type = "l",
       col = "blue",
       lwd = 2)
leg <- c("actual data", "smoothed estimate", "conditional sample")
legend(
  "topleft",
  leg,
  col = c("grey", "green", "blue"),
  lwd = c(1, 2, 2),
  cex = 0.7,
  bty = "n"
)
temp <- ts(c(out$epshat), start = 1871)
ylim <- range(out$epshat, eps.sim)
plot(
  temp,
  type = "l",
  col = "blue",
  ylim = c(-400, 300),
  ,
  main = "iii"
)
points(1871:1970,
       eps.sim,
       pch = 19,
       col = "sienna",
       cex = 0.5)
temp <- ts(c(out$etahat), start = 1871)
ylim <- range(out$etahat, eta.sim)
plot(
  temp,
  type = "l",
  col = "blue",
  ylim = c(-400, 300),
  ,
  main = "iv"
)
points(1871:1970,
       eta.sim,
       pch = 19,
       col = "sienna",
       cex = 0.5)

Nile.miss <- Nile
Nile.miss[21:40] <- NA
Nile.miss[61:80] <- NA
model.miss <-
  SSModel(Nile.miss ~ SSMtrend(1, Q = c(model$Q)) , H = c(model$H))
out.miss <- KFS(
  model.miss,
  filtering = 'state',
  smoothing = c('state', 'signal', 'disturbance'),
  simplify = F
)

par(mfrow = c(2, 2))
temp <- cbind(Nile.miss, c(out.miss$a)[-1])
cols <- c("grey", "blue")
plot.ts(
  temp,
  plot.type = "single",
  col = cols,
  xlab = "",
  ylab = "",
  main = "i",
  lwd = c(1, 2, 2)
)
legend(
  "topright",
  legend = c("Observation Data", "Filtered"),
  col = c("grey", "blue"),
  lwd = c(1, 2),
  bty = "n",
  cex = 0.8
)

temp <- ts(c(out.miss$P)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "ii"
)
temp <- cbind(Nile.miss, c(out.miss$alphahat))
cols <- c("grey", "blue")
plot.ts(
  temp,
  plot.type = "single",
  col = cols,
  xlab = "",
  ylab = "",
  main = "iii",
  lwd = c(1, 2, 2)
)
legend(
  "topright",
  legend = c("Observation Data", "Smoothed"),
  col = c("grey", "blue"),
  lwd = c(1, 2),
  bty = "n",
  cex = 0.8
)
temp <- ts(c(out.miss$V)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "iv"
)

#2.7 Forecasting for fig 2.6 -------

data(Nile)
Nile
Nile2 <- c(Nile, rep(NA, 30))

model <-
  SSModel(Nile2 ~ SSMtrend(1, Q = c(model$Q)) , H = c(model$H))
out <- KFS(
  model,
  filtering = 'state',
  smoothing = c('state', 'signal', 'disturbance'),
  simplify = F
)

par(mfrow = c(2, 2))
conf.bands <- cbind(out$a, as.vector(out$a) +
                      (sqrt(cbind(out$P)) * qnorm(0.95)) %*% (t(c(-1, 1))))
temp <- cbind(Nile2, conf.bands[-1,])
cols <- c("grey", "blue", "red", "red")
lwd <- c(1, 2, 1, 1)
lty <- c(1, 1, 2, 2)
plot.ts(
  temp,
  plot.type = "single",
  col = cols,
  lwd = lwd,
  lty = lty,
  xlab = "",
  ylab = "",
  main = "i"
)
legend(
  "topright",
  legend = c("Observation Data", "Forecasts", "Confidence intervals"),
  col = c("grey", "blue", "red"),
  lty = c(1, 1, 2),
  lwd = c(1, 2, 1),
  bty = "n",
  cex = 0.8
)
temp <- ts(c(out$P)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "ii"
)
temp <- ts(c(out$alphahat)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  ylim = c(700, 1250),
  main = "iii"
)
temp <- ts(c(out$F)[-1], start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "iv"
)
# 2.9 Param Est
model <-
  SSModel(Nile ~ SSMtrend(1, Q = list(matrix(NA))) , H = matrix(NA))
model <- fitSSM(inits = c(log(var(Nile)) , log(var(Nile))) ,
                model = model,
                method = 'BFGS')$model
c(model$H)
c(model$Q)

# 2.11 Diagnostic Checking for fig 2.7

out <-
  KFS(
    model,
    filtering = 'state',
    smoothing = c('state', 'signal', 'disturbance'),
    simplify = F
  )
stdres <- rstandard(out)
par(mfrow = c(2, 2))
temp <- ts(stdres, start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "i"
)
abline(h = 0, col = "grey")
par(mfrow = c(2, 2))
temp <- ts(stdres, start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "i"
)
abline(h = 0, col = "grey")

# 2.11 for fig 2.8
out <-
  KFS(
    model,
    filtering = 'state',
    smoothing = c('state', 'signal', 'disturbance'),
    simplify = F
  )
temp1 <- out$epshat / sqrt(c(out$V_eps))
temp2 <- out$etahat / sqrt(c(out$V_eta))

par(mfrow = c(2, 2))
temp <- ts(temp1, start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "i"
)
abline(h = 0, col = "grey")
hist(
  temp,
  prob = T,
  col = "grey",
  ylim = c(0, 0.3),
  main = "ii",
  xlab = "",
  ylab = ""
)
lines(density(temp), col = "blue", lwd = 2)

par(mfrow = c(2, 2))
temp <- ts(temp1, start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "",
  main = "i"
)
abline(h = 0, col = "grey")
hist(
  temp,
  prob = T,
  col = "grey",
  ylim = c(0, 0.3),
  main = "ii",
  xlab = "",
  ylab = ""
)
lines(density(temp), col = "blue", lwd = 2)
temp <- ts(temp2, start = 1871)
plot(
  temp,
  col = "blue",
  lwd = 2,
  xlab = "",
  ylab = "" ,
  main = "iii"
)
abline(h = 0, col = "grey")
hist(
  temp,
  prob = T,
  col = "grey",
  main = "iv",
  ylim = c(0, 1.2),
  xlab = "",
  ylab = ""
)
lines(density(temp), col = "blue", lwd = 2)

