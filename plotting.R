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
  ggtitle("Alpha")+xlab("Index")+ylab("Alpha")

plot_filtered <- ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line(aes(y = H_filtered))+ 
  ggtitle("Filtered")+xlab("Index")+ylab("Filtered")

plot_smoothed <-ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line(aes(y=H_smoothed))+ 
  ggtitle("Smoothed")+xlab("Index")+ylab("Smoothed")

plot_transformed <- ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line()+ 
  ggtitle("Transformed")+xlab("Index")+ylab("Transformed")

plot_demeaned
plot_transformed
plot_alpha
plot_filtered
plot_smoothed
plot_transformed