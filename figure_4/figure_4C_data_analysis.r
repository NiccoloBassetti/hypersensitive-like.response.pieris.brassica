############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "figure_4C_dataset.txt"
###
### Author R script:   Niccolo` Bassetti - niccolo.bassetti@protonmail.com 
###
### Principal Investigator:   Nina Fatouros - nina.fatouros@wur.nl


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
x=c("ggplot2", "agricolae", "RColorBrewer", "car", "multcomp", "lme4", "lsmeans", "tidyverse", "binom", "dplyr")

lapply(x, install.packages, character.only = TRUE)

# Set up working directory
getwd()                 # see current working directory
setwd("paste_here_your_path_to_file")	  # set working directory where input datasets are located. 



#### Statistical analysis ####
#____________________________#
packages=c("multcomp", "multcompView", "FSA", "rcompanion")
lapply(packages, library, character.only = TRUE)

mydata <- read.table("figure_4C_dataset.txt",header=TRUE,sep="")
data  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
str(data)

# Tranform class of categorical variables
data$treatment = factor(data$treatment, levels = unique(data$treatment)) # helps to plot x-axis variable in correct order
data$Response = data$treatment # helps to plot x-axis variable in correct order
str(data)

data


# Test for normality and homogeneity of variances
qqnorm(log2(data$relative_expression))
qqline(log2(data$relative_expression))
shapiro.test(log2(data$relative_expression))		# test for normality (small P = NOT NORMAL)
fligner.test(log2(data$relative_expression)~data$Response)	# test for equal variances (small P = NOT EQUAL VARIANCE)
      
# linear model
lm1 <- lm(log2(relative_expression) ~ Response, data)

summary(lm1)
anova(lm1)
AIC(lm1) # goodness-of-fit measure - smaller values are better
BIC(lm1) # goodness-of-fit measure - smaller values are better
coef(lm1) # coefficients of the model
confint(lm1) # confidence interval



# multiple comparison
pairs <- glht(lm1, linfct = mcp(Response = "Tukey"))
cld(pairs, level=0.01)

summary(pairs) 
confint(pairs)



#### Plot results ####
#____________________#

packages = c("ggplot2", "svglite","agricolae","ggtext")
lapply(packages, library, character.only = TRUE)


# Reformat dataframe for plotting purpose
data  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
str(data)
data$treatment = factor(data$treatment, levels = unique(data$treatment)) # helps to plot x-axis variable in correct order
data$Response = data$treatment # helps to plot x-axis variable in correct order
data

# rearrange factors for plotting
data$treatment = factor(data$treatment, levels=rev(levels(data$treatment)), ordered = TRUE)

head(data)
tail(data)

str(data)

# Prepare color palette
display.brewer.all(colorblindFriendly = T)
cbPalette=c("Paired") # choose palette colours

n=10
brewer.pal(n = n, name = cbPalette)   # visualize color names
display.brewer.pal(n = n, name = cbPalette) # inspect colors/


# Boxplot
par(mfrow=c(1,1))		# open graphic device (row, column) 

size=20

plot = ggplot(data) + # use geom_boxplot or geom_violin
  geom_boxplot(aes(x=Response, y=relative_expression
                   #,fill=treatment
                   ,color=treatment
  ),
  #colour="black", # change to black if you need to plot >1 Treatments
  outlier.colour=NA) +
  geom_point(mapping=aes(x=Response,y=relative_expression
                         ,fill=treatment
                         #,color=Genotype
  ), 
  colour="black", # black outline for points
  size=2,
  shape=19,
  alpha=0.5,
  position=position_jitterdodge(jitter.width=0.15,jitter.height=0,dodge.width=0.75,seed=NA)) +
  #scale_fill_manual(values=rep("white",length(levels(data$Response)))) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_fill_brewer(palette=cbPalette) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_colour_brewer(palette=cbPalette) + # set color for line/points - 
  scale_colour_manual(values = c("gray", "#E31A1C", "#1F78B4")) + # set color for line/points - 
  theme_bw() + # remove gray background 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  labs(x = "HR-like phenotype",
       y = "ethylene [ppm/mg FW]",
       size=size) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust =0.5, size=size, color="black"),
        axis.text.y=element_text(size=size, color="black"),
        axis.title=element_text(size=size, color="black"),
        axis.title.x=ggtext::element_markdown(margin = margin(t = size)),
        axis.title.y=ggtext::element_markdown(margin = margin(r = size))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  #theme(legend.position="none") + # change to set position, "none" remove legend
  theme(legend.position = c(0.20, 0.90)) + # remove this to override "none"+
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size=(size-5), face="plain")) + # change "face" value for style text
  #scale_x_discrete(limits=c("no_HR", "HR")) + # reorder scale x and legend
  scale_y_continuous(breaks=c(0,seq(100,300, by=100))
  )+
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
plot(plot)
ggsave("figure_4C.png", plot=plot, width = 6, height = 6, units="in") # CHANGE
dev.off()
