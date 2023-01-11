### add some intro material here ###
install.packages("bestNormalize")
install.packages("vctrs")
install.packages("bbmle")
install.packages("AICcmodavg")
install.packages("summarytools")

library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(datarium)

library(reshape2)
library(ggpmisc)
library(splus2R)

library(nlme)
library(agricolae)
library(multcomp)




# get data and format #

setwd("~/Git/Clam-OA-project")
data <-read.csv("clam_full_data.csv")
colnames(data)[1] <- gsub('^...','',colnames(data)[1])

# add column for treatment/species groups #

data$ts <- substring(data$If,1,3)
data$ts <- factor(data$ts,     
                  levels = c("L-C", "L-T", "M-C", "M-T"))

library(dplyr)

data$treatment <- recode_factor(data$treatment, C = "Ambient", T = "OA")

my_comparisons <- list( c("L-C", "L-T"), c("M-C", "M-T"), c("L-C", "M-C"), c("L-T", "M-T") )


### normalize Respirometry data ###

# data.resp <- data
# data.clean.resp <- na.omit(data.resp)

# calculate clam volume

data$length <- as.numeric(data$length)
data$width <- as.numeric(data$width)
data$depth <- as.numeric(data$depth)

data$volume <- ((data$length/2) * (data$width/2) * (data$depth/2)/1000) * (4/3) * pi

# calculate clam mass based on volume (equations derived from excel) and normalize O2 

data$mass <-  ifelse(data$species == "littleneck",
                                (data$volume * 0.0518) - 0.0229, 
                                (data$volume * 0.0617) - 0.0288)

data$umol.O2.normal <- data$umol.O2/data$mass

#data.split.resp.ts <- split(data.clean.resp, data.clean.resp$ts, drop=TRUE)

#data.split.resp.s  <- split(data.clean.resp, data.clean.resp$species, drop=TRUE)

####### Statistics #######

# factorize the various variable # 

data$species <- factor(data$species, levels = c("manila","littleneck"))
data$sex <- factor(data$sex, levels = c("f","m"))
data$spawned <- factor(data$spawned, levels = c("n","y"))
data$stage <- factor(data$stage, levels = c("5","4","3"))
data$mort <- factor(data$mort, levels = c("n","y"))

### respiration ANOVAs ###

# check normality of data #

test_resp <- data$umol.O2.normal

qqnorm(test_resp, main = "Q-Q Plot: untransformed") # check linearity
qqline(test_resp)
norm_test <- shapiro.test(test_resp) # p-value > 0.05 = good, don't need transformation
print(paste("shapiro test p-value, untransformed:", norm_test$p.value))

# Normalize response variable if normality test failed (spoiler, it did)

library(bestNormalize)

if(norm_test$p.value<0.05)     {
  normalized <- bestNormalize(test_resp, main = "Q-Q Plot: transformed")
  test_resp <- normalized$x.t # overwrite
  qqnorm(test_resp) # check linearity of transformed response
  qqline(test_resp)
  norm_test <- shapiro.test(test_resp) # p-value > 0.05 = good
  print(paste("shapiro test p-value, transformed:", norm_test$p.value))}
data$response <- test_resp # overwrite

# run ANOVA #

my_test <- aov(response ~ species * treatment * sex * spawned * stage, data = data)
my_test_summary <- summary(my_test)
summary(my_test)

# Compare model AIC scores (lowest score wins)
other <- aov(response ~ species * sex * stage, data = data)
other2 <- aov(response ~ species * treatment * spawned, data = data)
other3 <- aov(response ~ species * spawned, data = data)
other4 <- aov(response ~ treatment * spawned, data = data)
other5 <- aov(response ~ treatment * sex * stage, data = data)
other6 <- aov(response ~ species * treatment * sex * stage, data = data)

model.set <- list(my_test, other, other2, other3, other4, other5, other6)
model.names <- c("species:treatment:sex:spawned:stage", "species:sex:stage","species:spawned:treatment",
                 "species:spawned","treatment:spawned","treatment:sex:stage","species:treatment:sex:stage")

library(AICcmodavg)

aictab(model.set, modnames = model.names)

other2_summary <- summary(other2)
summary(other2)

other4_summary <- summary(other4)
summary(other4)

### ATPase ###

test_atp <- data$ATPase

qqnorm(test_atp, main = "Q-Q Plot: untransformed") # check linearity
qqline(test_atp)
norm_test_atp <- shapiro.test(test_atp) # p-value > 0.05 = good, don't need transformation
print(paste("shapiro test p-value, untransformed:", norm_test_atp$p.value))

# run ANOVA #

my_test_atp <- aov(ATPase ~ species * treatment * sex * spawned * stage, data = data)
my_test_atp_summary <- summary(my_test_atp)
summary(my_test_atp)

# Compare model AIC scores (lowest score wins)
other <- aov(ATPase ~ sex * stage, data = data)
other2 <- aov(ATPase ~ treatment * spawned * stage, data = data)
other3 <- aov(ATPase ~ species * spawned * stage, data = data)
other4 <- aov(ATPase ~ sex * spawned * stage, data = data)

model.set <- list(my_test_atp, other, other2, other3, other4)
model.names <- c("species:treatment:sex:spawned:stage", "sex:stage",
                 "stage:spawned:treatment","species:spawned:stage","sex:spawned:stage")

aictab(model.set, modnames = model.names)


other_summary <- summary(other)
summary(other)

other2_summary <- summary(other2)
summary(other2)

other3_summary <- summary(other3)
summary(other3)

other4_summary <- summary(other4)
summary(other4)

### spawned ###

# Calculate the chi-squared statistic and p-value for the test
data$sp_sex <- paste(data$species,data$sex)

# Identify rows with complete cases (i.e., no NA values) in the `cyl` column
complete_rows <- complete.cases(data$stage)

# Subset the data frame to include only rows with complete cases in the `cyl` column
data.naomit <- data[complete_rows,]

manila_f <- subset(data, sp_sex == "manila f")
manila_m <- subset(data, sp_sex == "manila m")
littleneck_f <- subset(data, sp_sex == "littleneck f")
littleneck_m <- subset(data, sp_sex == "littleneck m")

manila_f_st <- subset(data.naomit, sp_sex == "manila f")
manila_m_st <- subset(data.naomit, sp_sex == "manila m")
littleneck_f_st <- subset(data.naomit, sp_sex == "littleneck f")
littleneck_m_st <- subset(data.naomit, sp_sex == "littleneck m")

chisq.test(manila_f$treatment, manila_f$stage)
chisq.test(manila_m$treatment, manila_m$stage)
chisq.test(littleneck_f$treatment, littleneck_f$stage)
chisq.test(littleneck_m$treatment, littleneck_m$stage)

chisq.test(manila_f$treatment, manila_f$spawned)
chisq.test(manila_m$treatment, manila_m$spawned)
chisq.test(littleneck_f$treatment, littleneck_f$spawned)
chisq.test(littleneck_m$treatment, littleneck_m$spawned)

chisq.test(manila_f$treatment, manila_f$mort)
chisq.test(manila_m$treatment, manila_m$mort)
chisq.test(littleneck_f$treatment, littleneck_f$mort)
chisq.test(littleneck_m$treatment, littleneck_m$mort)

chisq.test(manila_f_st$stage, manila_f$spawned)
chisq.test(manila_m_st$stage, manila_m$spawned)
chisq.test(littleneck_f_st$stage, littleneck_f$spawned)
chisq.test(littleneck_m_st$stage, littleneck_m$spawned)

library(summarytools)


# fourth method:

manila_f %$%
  ctable(treatment, spawned,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

manila_m %$%
  ctable(treatment, spawned,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

littleneck_f %$%
  ctable(treatment, spawned,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

littleneck_m %$%
  ctable(treatment, spawned,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

## morts ##

manila_f %$%
  ctable(treatment, mort,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

manila_m %$%
  ctable(treatment, mort,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

littleneck_f %$%
  ctable(treatment, mort,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

littleneck_m %$%
  ctable(treatment, mort,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

## stage ##

manila_f_st %$%
  ctable(treatment, stage,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

manila_m_st %$%
  ctable(treatment, stage,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

littleneck_f_st %$%
  ctable(treatment, stage,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

littleneck_m_st %$%
  ctable(treatment, stage,
         prop = "r", chisq = TRUE, headings = FALSE
  ) %>%
  print(
    method = "render",
    style = "rmarkdown",
    footnote = NA
  )

###### graphics #####

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
