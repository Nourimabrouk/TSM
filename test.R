  transformed_df = transform_data(stockdata, returns)
  input_stocks = transformed_df[[1]]
  input_returns = transformed_df[[2]]
  initial_parameters = initialise_parameters_QML()
  
  state_space_parameters = initial_parameters[[1]]
  par_ini = initial_parameters[[2]]
  
  QML_params_returns <- optimize_parameters(input_returns, par_ini, state_space_parameters, TRUE) # (Print_output = TRUE)
  outputKalman_returns <- compute_kalmanfilter(input_returns, QML_params_returns, state_space_parameters)
  outputSmooth_returns <- compute_smoothed_state(input_returns, QML_params_returns, outputKalman_returns)
  
  QML_params_stocks <- optimize_parameters(input_stocks, par_ini, state_space_parameters, TRUE)
  outputKalman_stocks <- compute_kalmanfilter(input_stocks[,1], QML_params_stocks, state_space_parameters)
  outputSmooth_stocks <- compute_smoothed_state(input_stocks[,1],QML_params_stocks, outputKalman_stocks)

transform_data <- function(stockdata, returns){
  y = diff(log(stockdata$Close))
  x <- log((y - mean(y))^2)
  rv <- stockdata$RV[-1]
  
  input_matrix_stocks <- cbind(x, rv)
  input_returns <- returns$transformed

  stock_data <- cbind(x, stockdata$RV[-1])
  ret_trans <- returns$transformed
  return(list(stock_data, ret_trans))
}
initialise_parameters_QML <- function(){
  
  sig_eps <- (pi^2)/2 # Given in assignment
  mean_u <- -1.27 # Given in assignment
  
  par_ini <- c(0.1082, 0.991, -0.207, 0.0)
  
  Beta <- 0
  
  state_space_parameters <- data.frame(
    Q = 1,
    Z = 1,
    H = sig_eps,
    R = par_ini[1],
    T = par_ini[2],
    c = par_ini[3],
    d = mean_u,
    Beta = par_ini[4]
  )
  return(list(state_space_parameters, par_ini))
}
