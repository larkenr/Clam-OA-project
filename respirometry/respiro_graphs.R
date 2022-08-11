# Title: generate plots for diploid and triploid pacific oysters temperature experiment
# Author: Matthew George; mattgeorgephd@gmail.com
# Date: 05/2021

## clear
rm(list=ls())

## Grab the WD from the file location
library(rstudioapi)
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); getwd()

## Load R packages
library(readxl)
library(ggplot2)
library(stringr)
library(tidyverse)
library(Johnson)
library(agricolae)
library(nlme)
library(multcomp)

## Set ggplot theme
my_theme <- theme(line              = element_line(size=1.5),
                  rect              = element_rect(size=1.5),
                  text              = element_text(size=14,color="black"),
                  panel.background  = element_blank(),
                  panel.grid.major  = element_blank(), 
                  panel.grid.minor  = element_blank(),
                  axis.text.x       = element_text(size=16,color="black"),
                  axis.text.y       = element_text(size=16,color="black"),
                  axis.title.x      = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.title.y      = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.ticks.x      = element_line(color="black"),
                  axis.ticks.y      = element_line(color="black"),
                  # axis.line         = element_line(color = "black", size = 0.1),
                  panel.border      = element_rect(color = "black", fill=NA, size=1.5),
                  legend.key        = element_blank()) # removes background of legend bullets

### [4] Boxplot - SMR vs. time by ploidy - control treatment

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('output'); getwd()

trt_list <- read_excel("processed_summary_less.xlsx", sheet = "trt_list", col_names = TRUE)

MR_plot           <- read_excel("processed_summary_less.xlsx", sheet = "SMR", col_names = TRUE)
#MR_plot$species    <- factor(MR_plot$species, levels=c("L","M"),ordered=TRUE)
MR_plot$trt       <- factor(MR_plot$trt,levels=trt_list$trt_list,ordered=TRUE)

bp1 <- ggplot(MR_plot, aes(x=trt_list, y=SMR, group=as.factor(trt_list), fill=species)) +
  geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
               outlier.size=1, notch=FALSE) +
  scale_fill_manual(values=c("royalblue1", "orangered1")) +
  # geom_point() +
  scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 42)) +
  # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) + 
  theme(line              = element_line(size=0.8),
        rect              = element_rect(size=1),
        text              = element_text(size=14,color="black"),
        panel.background  = element_blank(),
        panel.grid.major  = element_blank(), 
        panel.grid.minor  = element_blank(),
        axis.text.x       = element_text(size=16,color="black"),
        axis.text.y       = element_text(size=16,color="black"),
        axis.title.x      = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y      = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.ticks.x      = element_line(color="black"),
        axis.ticks.y      = element_line(color="black"),
        # axis.line         = element_line(color = "black", size = 0.1),
        panel.border      = element_rect(color = "black", fill=NA, size=1.5),
        legend.key        = element_blank()) # removes background of legend bullets

bp1

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('plots'); getwd()

ggsave("BOXPLOT_SMR_timeseries.tiff",
       plot   = bp1,
       dpi    = 600,
       device = "jpeg",
       width  = 9,
       height = 6,
       units  = "in")