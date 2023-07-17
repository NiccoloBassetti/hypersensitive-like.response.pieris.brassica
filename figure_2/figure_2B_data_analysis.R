############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "figure_2B_dataset.txt"
###
### Author R script:   Dr. Niccolo` Bassetti - niccolo.bassetti@protonmail.com 
###
### Principal Investigator:   Dr. Nina Fatouros - nina.fatouros@wur.nl


# Clean R workspace
rm(list=ls())
gc()

# set up new R workspace
sessionInfo() #check packages available

.libPaths()   # see Path to R library
library()     # inspect packages available in Path


# Install R packages required for this analysis
update.packages()   # update all packages already installed

# install multiple packages at same time

library(lme4); library(lmerTest); library(pbkrtest); library(MASS)
packages=c("ggplot2", "agricolae", "RColorBrewer", "car", "multcomp", "lme4", "lsmeans", "tidyverse", "binom", "dplyr", "multcompView")

lapply(packages, install.packages, character.only = TRUE)

lapply(packages, library, character.only = TRUE)


# Set up working directory
getwd()                 # see current working directory

setwd()	  # set working directory where input datasets are located.




#### Prepare dataset ####
#_______________________#

mydata <- read.table(file = "figure_2B_dataset.txt",sep = "\t", header = TRUE) # !!! test sep character

str(mydata)
head(mydata)


# Reformat dataframe for plotting purpose
data_temp  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
data_temp$Treatment = factor(data_temp$Treatment, levels = unique(data_temp$Treatment)) # helps to plot x-axis variable in correct order
data_temp


# Filter data to keep only Eggs or Egg wash treatment
head(data_temp)
tail(data_temp)

data = data_temp

str(data)


#### Statistical analysis ####
#____________________________#

packages=c("multcomp", "multcompView", "FSA", "rcompanion")

lapply(packages, library, character.only = TRUE)


# Test for normality and homogeneity of variances
qqnorm(data$hr_score)
qqline(data$hr_score)
shapiro.test(data$hr_score)		# test for normality (small P = NOT NORMAL)
fligner.test(data$hr_score~data$Treatment)	# test for equal variances (small P = NOT EQUAL VARIANCE)


# Transform class variables
sapply(FUN = class, data)
data$Treatment = factor(data$Treatment, levels = unique(data$Treatment))
str(data)


# non-parametric test
kruskal.test(data$hr_score~data$Treatment)

# Data were plotted in Excel to generate figure 2B.