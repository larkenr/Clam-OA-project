library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(datarium)
library(dplyr)

setwd("~/Git/Clam-OA-project")
data <-read.csv("ATPase data sheet.csv")
names(data) <- c("sample","ID","treatment","A","B","A-B","Protein","activity","Mg")

data.clean <- na.omit(data)
data.clean$treat <- substring(data.clean$treatment,3,3)
data.clean$species <- substring(data.clean$treatment,1,1)

my_comparisons <- list( c("L-C", "L-T"), c("M-C", "M-T"), c("L-C", "M-C"), c("L-T", "M-T") )

# density all data
ggdensity(data.clean, x = "activity",
          fill = "treat",
          add = "median")


ggdensity(data.clean, x = "activity",
          fill = "species",
          add = "median")

ggdensity(data.clean, x = "activity",
          fill = "treatment",
          add = "median") +

 facet_wrap(treatment~.) +
 theme(strip.background = element_blank())

ggviolin(data.clean, x = "treatment", y = "activity",
         fill = "treatment") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

as.factor(data.clean$treatment)

levels(data.clean$treatment)

group_by(data.clean,treatment) %>%
  summarise(
    count = n(),
    mean = mean(activity, na.rm = TRUE),
    sd = sd(activity, na.rm = TRUE),
    median = median(activity, na.rm = TRUE),
    IQR = IQR(activity, na.rm = TRUE)
  )

ggboxplot(data.clean, x = "treatment", y = "activity", 
          color = "treatment",
          ylab = "activity", xlab = "Treatment")

kruskal.test(activity ~ treatment, data = data.clean)

pairwise.wilcox.test(data.clean$activity, data.clean$treatment,
                     p.adjust.method = "BH")
