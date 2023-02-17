library(seacarb)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)

setwd("~/Git/Clam-OA-project")

data <-read.csv("enviro_data.csv")
colnames(data)[1] <- gsub('^...','',colnames(data)[1])
data$group <- paste(data$treatment,data$rep_type)

data_clean <- na.omit(data)

mean_data <- data_clean %>%
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
mean_data$dic <- carb(flag = 8, mean_data$ph_with_correction, mean_data$alkalinity/1000000, T = mean_data$spectrometery_temperature, S = mean_data$salinity)$DIC * 1000000

mean_data$insitu_pH <- carb(flag = 15, mean_data$alkalinity/1000000, mean_data$dic/1000000, T = mean_data$insitu_temperature, S = mean_data$salinity)$pH

isp <- ggplot(mean_data, aes(x=day,y=insitu_pH,color = group)) +
  geom_line() + 
  xlab("")
isp
