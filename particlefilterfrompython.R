particle_filter <- function(stockdata){
  set.seed(1233)
  
  y <- diff(log(stockdata$Close))
  x <- log((y - mean(y))^2)
  
  N = 10;  n = length(y); 
  omega = -0.088; phi = 0.991; sigma_eta = 0.084
  
  a = rep(0, 100)
  H = matrix(data = 0, nrow = N, ncol = n)
  
  xi = omega/(1-phi)
  
  H[,1] = rnorm(N, 0,sqrt(sigma_eta/(1-phi^2)))
  
  for (t in 2:n) {
    H[,t] = rnorm(N, phi*H[,t-1], sqrt(sigma_eta))
    w_tilde = dnorm(0,sqrt(exp(xi)*exp(H[,t]))) #likeliehood
    w = w_tilde / sum(w_tilde)
    
    a[t] = sum(w * H[,t])
    
    H = sample_n(data.frame(H), N, replace = TRUE, weight = w)
  }
  return(a)
}

particle_filtered_stock <- particle_filter(stockdata)
H_filtered_stock

length(a)

