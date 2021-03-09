rm(list=ls())
# Define 
# y_t
# x_i
# particle sets
# alpha_t_i
# Where do we use omega?

n = 100; 
omega = -0.088; phi = 0.991; sigma_eta = 0.084

theta_0 = initialise_theta(n, sigma_eta, phi)
first_draw = draw_values(theta_0, sigma_eta, phi)

y_t = rnorm(100)
x_i = rnorm(100)

normalised_weights = compute_normalised_weights(first_draw, y_t, x_i)

att_ptt = compute_att_ptt(normalised_weights, first_draw)
resampled_theta = resample_stratified(att_ptt)

#perform_particlefilter_routine(n, sigma_eta, phi, y_t)

perform_particlefilter_routine <- function(n, sigma_eta, phi, y_t){
  # Implementation from 14.5.3 DK
  # With bootstrap filter
  # And pg284 resampling
  
  #Sigma eta, phi from previous optimizations
  theta = initialise_theta(n, sigma_eta, phi) #(theta_0)
  for (t in 1:n) {
    # Initialise values <- current, draw_values
    # Generate vector (step i of recursion)
    rnorm(n, phi*theta, sigma_eta^2)
    values = draw_values(n, sigma_eta, phi) # Vector of length n : theta_0 as random sample from normal unconditional distribution of theta
    normalised_weights = compute_normalised_weights(sigma, values, y_t)
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

compute_normalised_weights<-function(theta_t, y_t, x_i){
  
  sigma = sqrt(exp(x_i)) # answer to lucas discussion board
  weights = exp ( -log (2*pi*sigma^2 ) / 2  - theta_t / 2 - ( exp(- theta_t) * y_t^2)) / (2 * sigma ^ 2) # 322 DK (ii)   
  normalised_weights = weights / sum ( weights )
  
  return(normalised_weights)
}
compute_att_ptt<- function(weights, theta){
  a_hat_t_t = sum ( weights * theta)
  p_hat_t_t = sum ( weights * theta ^ 2 - a_hat_t_t ^ 2 )
  
  return(data.frame(a_hat_t_t, p_hat_t_t, weights))
}

resample_stratified<- function(df_attptt){

  # Resampling from pg 284
    alpha_t = sample_n(df_attptt[,1:2], 100, replace = TRUE, weight = df_attptt[,3])
  
  return(alpha_t)
}
