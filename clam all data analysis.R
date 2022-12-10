install.packages("reshape2")
install.packages("ggpmisc")

library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(datarium)
library(dplyr)
library(reshape2)
library(ggpmisc)
library(splus2R)

setwd("~/Git/Clam-OA-project")
data <-read.csv("clam full data.csv")
data$ts <- substring(data$ID,1,3)
data$ts <- factor(data$ts,     
                  levels = c("L-C", "L-T", "M-C", "M-T"))
data$treatment <- substring(data$ID,3,3)
data$treatment <- recode_factor(data$treatment, C = "Ambient", T = "OA")
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

data.split.resp.ts <- split(data.clean.resp, data.clean.resp$ts, drop=TRUE)

data.split.resp.s  <- split(data.clean.resp, data.clean.resp$species, drop=TRUE)

# density graphs

ggdensity(data.split.resp.ts[["L-C"]], x = "umol.O2.normal",
#          fill = "sex",
          add = "median")

ggdensity(data.clean.resp, x = "umol.O2.normal",
          fill = "treatment",
          add = "median")

ggdensity(data.clean.resp, x = "umol.O2.normal",
          fill = "treatment") +
  xlab("Oxygen consumption (umol O2/L/hr/g)")    +
  facet_wrap(treatment ~ species) +
#  facet_wrap(ts~.) 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        text = element_text(size = 20),
        legend.title = element_blank())

data.clean.resp <- na.omit(data.clean.resp) 

data.clean.resp %>%
  group_by(treatment) %>%
  summarise(median(umol.O2.normal))
data.clean.resp$treatment <- as.character(data.clean.resp$treatment)

dl <- density(data.split.resp.ts[["L-C"]]$umol.O2.normal)

dm <- density(data.split.resp.s[["M"]]$umol.O2.normal)

localMaxima <- function(x) {
  y <- diff(c(-Inf,x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

localMinima <- function(x) {
  y <- diff(c(Inf,x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

loc.max.l <- dl$x[localMaxima(dl$y)]
loc.max.m <- dm$x[localMaxima(dm$y)]
loc.min.l <- dl$x[localMinima(dl$y)]
loc.min.m <- dm$x[localMinima(dm$y)]

sum(data.split.resp$`M-C`$umol.O2.normal < 4.39)
sum(data.split.resp$`M-C`$umol.O2.normal > 4.39)

ts <- c("L-C","L-T","M-C","M-T")
sleepy <- c(15,19,7,15)
active <- c(22,19,35,23)

activity.df <- data.frame(ts,sleepy,active)

ggbarplot(activity.df,x="ts", y=c("sleepy","active"),merge = TRUE)

# violin and boxplots

ggviolin(data.clean.resp, x = "ts", y = "umol.O2.normal",
         fill = "ts") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")


ggboxplot(data.clean.resp, x = "stage", y = "umol.O2.normal",
          fill = "treatment") +
  facet_wrap(ts~.)            +
  ylab("Oxygen consumption (umol O2/L/hr/gr)")    +
  xlab("")                    +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 7) +
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        axis.text.x = element_blank())

#  stats

kruskal.test(umol.O2.normal ~ ts, data = data.clean.resp)

pairwise.wilcox.test(data.clean.resp$umol.O2.normal, data.clean.resp$spawned,
                     p.adjust.method = "BH")


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
          fill = "ts",
          add = "median")  +
  xlab("ATPase activity")  +
  facet_wrap(ts~.) +
  theme(strip.background = element_blank())

ggviolin(data.clean.atp, x = "ts", y = "ATPase",
         fill = "treatment") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

ggboxplot(data.clean.atp, x = "sex", y = "ATPase",
          fill = "treatment") +
  facet_wrap(ts~.)            +
  ylab("ATPase activity")     +
  xlab("")                    +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 7) +
  theme(text = element_text(size = 20),
        legend.title = element_blank())
#        axis.text.x = element_blank())

kruskal.test(ATPase ~ stage, data = data.clean.atp)

pairwise.wilcox.test(data.clean.atp$ATPase, data.clean.atp$stage,
                     p.adjust.method = "BH")

### ATPase and Respiration

ggscatter(data.clean.resp,x="ATPase",y="stage",color="species",palette=c("red","black"))


### 

data.stage <- data
data.clean.stage <- na.omit(data.stage)
data.clean.stage$stage <- as.character(data.clean.stage$stage)
data.clean.stage$sts <- paste(data.clean.stage$ts,data.clean.stage$sex,sep="-")

###



ggplot(data.clean.stage, aes(x = sts, fill = stage)) +
  geom_bar(position = "fill", stat = "count") +
  geom_text(aes(label = paste0("n=", after_stat(count))), stat='count', position = position_fill(vjust = 0.5)) +
  theme_classic()
