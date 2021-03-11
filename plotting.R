#rm(list=ls())

library(ggplot2)
library(ggthemes)

plot_returns_input = readRDS(file = "plot_returns_input.rds")
plot_stock_input = readRDS(file = "plot_stock_input.rds")


plot_demeaned <- ggplot(plot_returns_input, aes(index, demeaned))+
  theme_minimal()+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_y_continuous(breaks = seq(-.05,.05,.025))+
  scale_x_continuous(breaks = seq(0,900,100))+
  ggtitle("Demeaned")+xlab("Index")+ylab("Demeaned")

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


plot_d_combined <-ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line(aes(y=H_smoothed))+
  geom_line(aes(y=H_filtered), col = "red") +
  ggtitle("Filtered Ht and Smoothed Ht")+xlab("Index")+ylab("Smoothed")


  ggtitle("Smoothed")+xlab("Index")+ylab("Smoothed")

plot_d_combined <-ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line(aes(y=H_smoothed))+
  geom_line(aes(y=H_filtered), col = "red") +
  ggtitle("Smoothed")+xlab("Index")+ylab("Smoothed")

plot_d_combined
plot_demeaned
plot_transformed
plot_alpha
plot_filtered
plot_smoothed
plot_d_combined

#plot 2e
colnames(plot_returns_input)
colnames(plot_stock_input)

# 
plot_one <-ggplot(plot_stock_input, aes(x = index))+
theme_minimal()+
  geom_line(aes(y=alpha), col = "blue")+
  geom_line(aes(y=alpha_rv), col = "red") +
  geom_point(aes(y = x), size = 0.1)+
  ggtitle("Plot 1 zeus")+xlab("Index")+ylab("")

plot_two <-ggplot(plot_stock_input, aes(x = index))+
theme_minimal()+
  geom_line(aes(y=H_filtered_stock))+
  geom_line(aes(y=H_filtered_stock_rv), col = "red") +
  ggtitle("Plot 2 zeus")+xlab("Index")+ylab("")

plot_three <- ggplot(plot_stock_input, aes(x = index))+
  theme_minimal()+
    geom_line(aes(y=H_smoothed))+
    geom_line(aes(y=H_smoothed_rv), col = "red") +
    ggtitle("Plot 3 zeus")+xlab("Index")+ylab("")

plot_particle <- ggplot(plot_stock_input, aes(x = index))+
  theme_minimal()+
  geom_line(aes(y=H_filtered_stock))+
  geom_line(aes(y=particle_filtered_stock), col = "red") +
  ggtitle("Particle filter plot")+xlab("Index")+ylab("Filtered / Particle filtered (red)")

plot_particle
