install.packages("ggpubr")
install.packages("rstatix")
install.packages("datarium")



library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(datarium)

data <-read.csv("respirometry no constraints.csv")
names(data) <- c("sample","umol","length","width","depth","volume")

data.clean <- na.omit(data)
data.clean$ts <- substring(data.clean$sample,1,3)
data.clean$treatment <- substring(data.clean$sample,3,3)
data.clean$species <- substring(data.clean$sample,1,1)
data.clean$umol.vol <- data.clean$umol/data.clean$volume
data.clean$scale.umol <- scale(data.clean$umol)
data.clean$scale.umol.vol <- scale(data.clean$umol.vol)
data.clean.cut <- subset(data.clean, data.clean$umol > 0 & data.clean$umol < 40)
cut.below <- subset(data.clean, data.clean$umol < 5)
cut.above <- subset(data.clean, data.clean$umol > 40)
cut <- rbind(cut.above, cut.below)
data.clean.cut$scale.umol <- scale(data.clean.cut$umol)
data.clean.cut$scale.umol.vol <- scale(data.clean.cut$umol.vol)

my_comparisons <- list( c("L-C", "L-T"), c("M-C", "M-T"), c("L-C", "M-C"), c("L-T", "M-T") )

# all data
ggdensity(data.clean, x = "umol",
          fill = "treatment",
          add = "median")


ggdensity(data.clean, x = "umol",
          fill = "species",
          add = "median")

ggdensity(data.clean, x = "umol",
          fill = "ts",
          add = "median") +
  facet_wrap(ts~.) +
  theme(strip.background = element_blank())

ggdensity(data.clean, x = "scale.umol",
          fill = "ts",
          add = "median")



# cut data
ggdensity(data.clean.cut, x = "umol",
          fill = "treatment",
          add = "median")

ggdensity(data.clean.cut, x = "umol",
          fill = "species",
          add = "median")

ggdensity(data.clean.cut, x = "umol",
          fill = "ts",
          add = "median") +
  facet_wrap(.~ts) +
  theme(strip.background = element_blank())

ggdensity(data.clean.cut, x = "scale.umol",
          fill = "ts",
          add = "median")

ggviolin(data.clean.cut, x = "ts", y = "umol",
          fill = "ts") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")


# corrected by vol

ggdensity(data.clean, x = "scale.umol.vol",
          fill = "treatment",
          add = "median")

ggdensity(data.clean, x = "scale.umol.vol",
          fill = "species",
          add = "median")

ggdensity(data.clean, x = "scale.umol.vol",
          fill = "ts",
          add = "median") +
  facet_wrap(.~ts) +
  theme(strip.background = element_blank())

ggdensity(data.clean, x = "scale.umol.vol",
          fill = "ts",
          add = "median")

ggviolin(data.clean, x = "ts", y = "scale.umol.vol",
         fill = "ts") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

#cut and corrected by volume

ggdensity(data.clean.cut, x = "scale.umol.vol",
          fill = "treatment",
          add = "median")

ggdensity(data.clean.cut, x = "scale.umol.vol",
          fill = "species",
          add = "median")

ggdensity(data.clean.cut, x = "scale.umol.vol",
          fill = "ts",
          add = "median") +
  facet_wrap(.~ts) +
  theme(strip.background = element_blank())

ggdensity(data.clean.cut, x = "scale.umol.vol",
          fill = "ts",
          add = "median")

ggviolin(data.clean.cut, x = "ts", y = "scale.umol.vol",
         fill = "ts") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")
