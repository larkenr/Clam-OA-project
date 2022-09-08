install.packages("reshape2")

library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(datarium)
library(dplyr)
library(reshape2)

setwd("~/Git/Clam-OA-project")
data <-read.csv("clam full data.csv")
data$ts <- substring(data$ID,1,3)
data$ts <- factor(data$ts,     
                             levels = c("L-C", "L-T", "M-C", "M-T"))
data$treatment <- substring(data$ID,3,3)
data$species <- substring(data$ID,1,1)

my_comparisons <- list( c("L-C", "L-T"), c("M-C", "M-T"), c("L-C", "M-C"), c("L-T", "M-T") )

### Respirometry

data.resp <- data
data.clean.resp <- na.omit(data.resp)

# calculate clam volume

data.clean.resp$length <- as.numeric(data.clean.resp$length)
data.clean.resp$width <- as.numeric(data.clean.resp$width)
data.clean.resp$depth <- as.numeric(data.clean.resp$depth)

data.clean.resp$volume <- ((data.clean.resp$length/2) * (data.clean.resp$width/2) * (data.clean.resp$depth/2)/1000) * (4/3) * pi

# calculate clam mass based on volume (equations derived from excel) and normalize O2 

data.clean.resp$mass <-  ifelse(data.clean.resp$species == "L",
                                (data.clean.resp$volume * 0.0518) - 0.0229, 
                                (data.clean.resp$volume * 0.0617) - 0.0288)

data.clean.resp$umol.O2.normal <- data.clean.resp$umol.O2/data.clean.resp$mass

# density graphs

ggdensity(data.clean.resp, x = "umol.O2.normal",
          fill = "sex",
          add = "median")

ggdensity(data.clean.resp, x = "umol.O2.normal",
          fill = "species",
          add = "median")

ggdensity(data.clean.resp, x = "umol.O2.normal",
          fill = "treatment",
          add = "median")  +
  xlab("umol 02 per gram") +
  facet_wrap(ts~.) +
  theme(strip.background = element_blank())

ggviolin(data.clean.resp, x = "ts", y = "umol.O2.normal",
         fill = "ts") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

kruskal.test(umol.O2.normal ~ ts, data = data.clean.resp)

pairwise.wilcox.test(data.clean.resp$umol.O2.normal, data.clean.resp$ts,
                     p.adjust.method = "BH")

ggboxplot(data.clean.resp, x = "ts", y = "umol.O2.normal",
          fill = "treatment") +
  ylab("umol O2 per gram")    +
  xlab("")                    +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

### ATPase activity

data.atp <- data
data.atp$ATPase <- as.numeric(data.atp$ATPase)
data.clean.atp <- na.omit(data.atp)

# density graphs

ggdensity(data.clean.atp, x = "ATPase",
          fill = "sex",
          add = "median")

ggdensity(data.clean.atp, x = "ATPase",
          fill = "species",
          add = "median")

ggdensity(data.clean.atp, x = "ATPase",
          fill = "treatment",
          add = "median")  +
  xlab("ATPase activity")  +
  facet_wrap(ts~.) +
  theme(strip.background = element_blank())

ggboxplot(data.clean.atp, x = "ts", y = "ATPase",
          fill = "treatment") +
  ylab("ATPase activity")    +
  xlab("")                    +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

kruskal.test(ATPase ~ sex, data = data.clean.atp)

pairwise.wilcox.test(data.clean.atp$ATPase, data.clean.atp$ts,
                     p.adjust.method = "BH")

### Reproductive state


### 