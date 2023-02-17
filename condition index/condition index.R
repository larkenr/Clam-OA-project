library(ggplot2)

setwd("~/Git/Clam-OA-project/condition index")

data <-read.csv("condition_index.csv")

data$condition <- (data$dry/data$shell)*100

ggboxplot(data,"date","condition")

