

perform_particlefilter_routine <- function(n, sigma_eta, phi, sigma, theta_t, y_t){
  # Implementation from 14.5.3 DK
  # With bootstrap filter
  # And pg284 resampling

  
for (t in 1:n) {
  values = draw_values(n, sigma_eta, phi) # Vector of length n : theta_0 as random sample from normal unconditional distribution of theta
  normalised_weights = compute_normalised_weights(sigma, theta_t, y_t)
  c(att, ptt) = compute_att_ptt(theta_t, normalised_weights)
  a_t = resampling(att)
  return(a_t)
  } 
}

draw_values <- function(n, sigma_eta, phi){
  var = sigma_eta^2 / (1 - phi^2)
  theta_0 = rnorm(n, 0, var)
  return(theta_0)
}
compute_normalised_weights<- function(){
  weights = exp ( -log (2*pi*sigma^2 ) / 2  - theta_t / 2 - ( exp(- theta_t) * y_t^2)) / (2 * sigma ^ 2) # 322 DK (ii)   
  normalised_weights = weights / sum ( weights )
  return(normalised_weights)
}

compute_att_ptt<- function(weights, theta){
  a_hat_t_t = sum ( weights * theta_t_i)
  p_hat_t_t = sum ( weights * theta_t_i ^ 2 - a_hat_t_t ^ 2 )
  
  return(cbind(a_hat_t_t, p_hat_t_t, weights))
}

# placeholders
a_hat_t_t = 1:100
p_hat_t_t = 1:100
weights = 100:1 / sum(100:1)
df_att_ptt_weights <- cbind(a_hat_t_t, p_hat_t_t, weights)

resampling<- function(df_att_ptt_weights){
  # Resampling from pg 284

  att = df_att_ptt_weights[,1]
  weights = df_att_ptt_weights[,3]
  weights = abs(rnorm(100,1,.5))
  resampled = sample_n(input, 100, replace = TRUE, weight = input$weights)
  
  return(resampled)
}
