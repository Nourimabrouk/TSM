# Data prep
rm(list=ls())

library(ggplot2)
library(ggthemes)

plot_returns_input = readRDS(file = "plot_returns_input.rds")
plot_returns_input

# #Quickplots
# autoplot(plot_returns_input, demeaned)
# autoplot(plot_returns_input, transformed)
# autoplot(stonkdata, Close)
# autoplot(stonkdata, RV)
# range(plot_returns_input$demeaned)

# 14.5 (i) # Demeaned returns
plot_demeaned <- ggplot(plot_returns_input, aes(index, demeaned))+
  theme_minimal()+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_y_continuous(breaks = seq(-.05,.05,.025))+
  scale_x_continuous(breaks = seq(0,900,100))+ 
  ggtitle("Demeaned")+xlab("Index")+ylab("Demeaned")
# Transformed
plot_transformed <- ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line()+ 
  ggtitle("Transformed")+xlab("Index")+ylab("Transformed")

plot_alpha <- ggplot(data = plot_returns_input, aes(x = index)) +
  theme_minimal()+
  geom_line(aes(y = alpha), size = 2, color = "red")+
  geom_point(aes(y = transformed))+
  scale_x_continuous(breaks = seq(0,900,100)) + 
  ylim(-20,-5)+
  ggtitle("Smoothed value of h_t")+xlab("Index")+ylab("h_t")

# plot_filtered <- ggplot(plot_returns_input, aes(index, transformed))+
#   theme_minimal()+
#   geom_line(aes(y = H_filtered))+ 
#   ggtitle("Filtered Ht")+xlab("Index")+ylab("Filtered")
# 
# plot_smoothed <-ggplot(plot_returns_input, aes(index, transformed))+
#   theme_minimal()+
#   geom_line(aes(y=H_smoothed))+ 
#   ggtitle("Smoothed Ht")+xlab("Index")+ylab("Smoothed")

plot_d_combined <-ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line(aes(y=H_smoothed))+
<<<<<<< Updated upstream
  geom_line(aes(y=H_filtered), col = "red") + 
  ggtitle("Filtered Ht and Smoothed Ht")+xlab("Index")+ylab("Smoothed")


=======
  ggtitle("Smoothed")+xlab("Index")+ylab("Smoothed")

plot_d_combined <-ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line(aes(y=H_smoothed))+
  geom_line(aes(y=H_filtered), col = "red") + 
  ggtitle("Smoothed")+xlab("Index")+ylab("Smoothed")
  
plot_d_combined
>>>>>>> Stashed changes
plot_demeaned
plot_transformed
plot_alpha
plot_filtered
plot_smoothed
<<<<<<< Updated upstream
plot_d_combined
=======
>>>>>>> Stashed changes

#plot 2e

h_t_stock <- outputKalman_stocks$h_t
h_t_stock_rv <- outputKalman_stocks_rv$h_t

H_filtered_stock <- h_t_stock - xi_stocks
H_filtered_stock_rv <- h_t_stock_rv - xi_stocks_rv

plot(ts(H_filtered_stock), col="red", plot.type="single", ylab="", main="H_t Filtered")
lines(ts(H_filtered_stock_rv))

#plot 3e
H_smoothed <- outputSmooth_stocks$alpha - xi_stocks
H_smoothed_rv <- outputSmooth_stocks_rv$alpha - xi_stocks_rv

plot(ts(H_smoothed), col="red", plot.type="single", ylab="", main="H_t Smoothed")
plot(ts(H_smoothed_rv), col="red", plot.type="single", ylab="", main="H_t Smoothed")


plot(ts(outputSmooth_returns$alpha) , col="red", plot.type="single", ylab="", main="h_t", ylim=c(min(returns$transformed), max(returns$transformed)))
points(returns$transformed, col="black")

