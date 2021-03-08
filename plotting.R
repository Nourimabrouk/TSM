# Data prep

#Quickplots
autoplot(returns, demeaned)
autoplot(returns, transformed)

autoplot(stonkdata, Close)
autoplot(stonkdata, RV)
range(returns$demeaned)

# 14.5 (i)
ggplot(returns, aes(index, demeaned))+
  theme_minimal()+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_y_continuous(breaks = seq(-.05,.05,.025))+
  scale_x_continuous(breaks = seq(0,900,100))

# # 14.5 (ii)
# ggplot(returns, aes(index))+
#   theme_minimal()+
#   geom_point()+
#   # + geom_line(smoothed_estimate k+theta)
#   geom_hline(yintercept = 0)+
#   scale_x_continuous(breaks = seq(0,900,100))

# 14.5 (iii)
# ggplot(returns, aes(index, SE_volmeasure))+
#   theme_minimal()+
#   geom_line()


plotReturns<- function(svData, label){
  
}

plot2Lines<- function(kalman, kalmanLabel,data, dataLabel){
  
}

plotForecast<- function(DF_forecast, DF_Kalman, NileData_array, yearIncludingForecast, T){
  
}
plotKalmanSV<- function(){
  
}
