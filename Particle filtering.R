perform_particlefilter_routine <- function(input){
  # Implementation from 14.5.3 DK
  # With bootstrap filter
  # And ... resampling

  n = 100
  sigma_eta = .5
  phi = .5
  
for (t in 1:n) {
  values = draw_values(n, sigma_eta, phi) # Vector of length n : theta_0 as random sample from normal unconditional distribution of theta
  normalised_weights = compute__normalised_weights(sigma, theta_t, y_t)
  c(att, ptt) = compute_att_ptt(theta_t, normalised_weights)
  a_t_i = resampling()
  } 
}

draw_values <- function(n, sigma_eta, phi){
  var = sigma_eta^2 / (1 - phi^2)
  theta_0 = rnorm(n, 0, var)
  return(theta_0)
}
compute__normalised_weights<- function(){
  weights = exp ( -log (2*pi*sigma^2 ) / 2  - theta_t / 2 - ( exp(- theta_t) * y_t^2)) / (2 * sigma ^ 2) # 322 DK (ii)   
  normalised_weights = weights / sum ( weights )
  return(normalised_weights)
}

compute_att_ptt<- function(weights, theta){
  a_hat_t_t = sum ( weights * theta_t_i)
  p_hat_t_t = sum ( weights * theta_t_i ^ 2 - a_hat_t_t ^ 2 )
  
  return(c(a_hat_t_t, p_hat_t_t))
}
resampling<- function(list_att_ptt){
  # stratified sampling 12.3.3 ? 
  
  # Resampling from pg 281 or 284
  
  # resample with replacements and weights
  input = data.frame(
    sample = 1:100)
  weights = abs(rnorm(100,1,.5))
  resampled = sample_n(input, 100, replace = TRUE, weight = input$weights)

  resampled
}