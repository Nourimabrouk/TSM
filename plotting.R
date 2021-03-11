#rm(list=ls())
diff_rvbase = mean(plot_stock_input$H_filtered_stock - plot_stock_input$H_filtered_stock_rv)

diff_three = mean(plot_stock_input$H_smoothed - plot_stock_input$H_smoothed_rv)

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
  xlab("Time index")+ylab("Demeaned returns")
plot_transformed <- ggplot(plot_returns_input, aes(index, transformed))+
  theme_minimal()+
  geom_line()+
  xlab("Time index")+ylab("Transformed returns")
plot_alpha <- ggplot(data = plot_returns_input, aes(x = index)) +
  theme_minimal()+
  geom_line(aes(y = alpha, color = "alpha"), size = .85)+ # Add legend
  geom_point(aes(y = transformed, color = "transformed"), size = .75)+ # Add legend
  labs(x = "Time index",
       y = "",
       color = "Legend") +
  scale_x_continuous(breaks = seq(0,900,100)) +
  ylim(-20,-5)+xlab("Time index")+ylab("")

plot_d_combined <-ggplot(plot_returns_input, aes(index, transformed)) +
  # Drop first obs
  # Add legend smoothed / filtered
  theme_minimal() +
  geom_line(aes(y=H_smoothed)) + # add legend
  geom_line(aes(y=H_filtered), col = "red") + # add legend
  xlab("Index") + ylab("")

plot_one <- ggplot(plot_stock_input %>% slice(-c(1:3)), aes(x = index))+
  # Add legend
  # SPX data
  # Dots = Log x^2
  # Blue = smoothed ht with rv
  # Red = smoothed ht base model
  theme_minimal()+
  geom_line(aes(y=alpha), col = "blue", size=0.55)+
  geom_line(aes(y=alpha_rv), col = "red", size=0.55) +
  geom_point(aes(y = x), size = 0.1)+
  ggtitle("SPX log(x^2) with smoothed estimates h_t")+xlab("Time")+ylab("")

plot_two <-ggplot(plot_stock_input %>% slice(-c(1:2)), aes(x = index))+
theme_minimal()+
  geom_line(aes(y=H_filtered_stock), col = "blue")+
  geom_line(aes(y=H_filtered_stock_rv+diff_rvbase), col = "red") +
 xlab("Time index")+ylab("")

plot_three <- ggplot(plot_stock_input %>% slice(-c(1:3)), aes(x = index))+
    theme_minimal()+
    geom_line(aes(y=H_smoothed), col = "blue") +
    geom_line(aes(y=H_smoothed_rv + diff_three), col = "red") +
    xlab("Time index")+ylab("")

plot_particle <- ggplot(plot_stock_input, aes(x = index))+
  theme_minimal()+
  geom_line(aes(y=H_filtered_stock))+
  geom_line(aes(y=particle_filtered_stock), col = "red") +
  xlab("Time index")+ylab("Filtered / Particle filtered (red)")

plot_demeaned
plot_transformed
plot_alpha
plot_d_combined
plot_one
plot_two
plot_three
plot_particle