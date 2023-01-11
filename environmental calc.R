library(seacarb)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)

setwd("~/Git/Clam-OA-project")

data <-read.csv("enviro_data.csv")
colnames(data)[1] <- gsub('^...','',colnames(data)[1])
data$group <- paste(data$treatment,data$rep_type)

mean_data <- data %>%
  group_by(day,group) %>%
  summarize_all(mean)

### graph main parameters ###

#data <- data.frame(
#  day = as.Date("2017-06-14") - 0:364,
#  value = runif(365) + seq(-140, 224)^2 / 10000
#)

# Most basic bubble plot
s <- ggplot(mean_data, aes(x=day, y=salinity, color = group)) +
  geom_line() + 
  xlab("")
s

p <- ggplot(mean_data, aes(x=day,y=ph_with_correction,color = group)) +
  geom_line() + 
  xlab("")
p

t <- ggplot(mean_data, aes(x=day, y=insitu_temperature, color =group)) +
  geom_line() + 
  xlab("")
t

plot_grid(s,p,t,
          nrow=3)

### alkalinity calc ###
data$dic <- carb(flag = 8, data$ph_with_correction, data$alkalinity/1000000, T = data$spectrometery_temperature, S = data$salinity)$DIC * 1000000

data$insitu_pH <- carb(flag = 15, data$alkalinity/1000000, data$dic/1000000, T = data$insitu_temperature, S = data$salinity)$pH
