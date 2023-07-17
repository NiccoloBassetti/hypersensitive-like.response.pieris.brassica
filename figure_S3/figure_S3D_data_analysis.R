############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "figure_S2D_dataset.txt"
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
x=c("ggplot2", "agricolae", "RColorBrewer", "car", "multcomp", "lme4", "lsmeans", "tidyverse", "binom", "dplyr")

lapply(x, install.packages, character.only = TRUE)

# Set up working directory
getwd()                 # see current working directory
setwd("paste_here_your_path_to_file")	  # set working directory where input datasets are located. 



#### Prepare dataset ####
#_______________________#
mydata <- read.table("figure_S2D_dataset.txt",header=TRUE,sep="")

head(mydata)
names(mydata)

dim(mydata)
str(mydata)

# change class of dependent variables
mydata$Time_point <- as.factor(data$Timepoint)
mydata$Treatment <- as.factor(data$Treatment)



#### Statistical analysis ####
#____________________________#

packages=c("multcomp", "multcompView", "FSA", "rcompanion")
lapply(packages, library, character.only = TRUE)

# Check normality of the mydata
qqnorm(mydata$Relative_exp)
qqline(mydata$Relative_exp)
shapiro.test(mydata$Relative_exp)		# test for normality (small P = NOT NORMAL)
fligner.test(mydata$Relative_exp~mydata$Timepoint)	# test for equal variances (small P = NOT EQUAL VARIANCE)

# CONCLUSION
# Data non normal, a non-parametric test was used


# non-parametric test
kruskal.test(mydata$Relative_exp~mydata$Timepoint)

# pairwise comparison 1: Dunn test
PT = dunnTest(mydata$Relative_exp~mydata$Timepoint,
              data=mydata,
              method="bh")    # require "FSA", p.adjust to be used to adjust p-values
PT
PT = PT$res

cldList(comparison = PT$Comparison, # gives weird results, not correspondin with Wilcoxon
        p.value    = PT$P.unadj,
        threshold  = 0.05) # require "rcompanion"


#### Plot results ####
#____________________#

packages = c("ggplot2", "svglite","agricolae","ggtext")
lapply(packages, library, character.only = TRUE)

# Reformat dataframe for plotting purpose
data  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
data$Timepoint = factor(data$Timepoint, levels = unique(data$Timepoint)) # helps to plot x-axis variable in correct order

head(data)


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
  geom_boxplot(aes(x=Timepoint, y=Relative_expression
                   ,fill=Genotype
                   #,color="Genotype"
                  ),
                colour="#E31A1C", # change to black if you need to plot >1 Treatments
                outlier.colour=NA) +
  geom_point(mapping=aes(x=Timepoint,y=Relative_expression
                         ,fill=Genotype
                         #,color=Genotype
                         ), 
              colour="black", # black outline for points
              size=2,
              shape=19,
              alpha=0.5,
              position=position_jitterdodge(jitter.width=0.15,jitter.height=0,dodge.width=0.75,seed=NA)) +
  labs(x="", y="Relative  expression", size=size) +
  scale_fill_manual(values=rep("white",length(levels(data$Timepoint)))) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_fill_brewer(palette=cbPalette) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_colour_brewer(palette=cbPalette) + # set color for line/points - 
  theme_bw() + # remove gray background 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  labs(x = "Time (hours)",
       y = "*PR1* relative expression",
      ) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust =0.5, size=size, color="black"),
        axis.text.y=element_text(size=size, color="black"),
        axis.title=element_text(size=size, color="black"),
        axis.title.x=element_text(margin = margin(t = size)),
        axis.title.y = ggtext::element_markdown(margin = margin(r = size))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  theme(legend.position="none") + # change to set position, "none" remove legend
  #theme(legend.position = c(0.85, 0.90)) + # remove this to override "none"
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size=(size), face="plain")) + # change "face" value for style text
  #scale_x_discrete(limits=c("0h","3h","6h","24h","48h")) + # reorder scale x and legend
  scale_y_continuous(#breaks=c(0,seq(100,300, by=100)),
                      breaks=c(seq(0,12.5, by=2.5)),
                      limits=c(0,12.5)
                      )+
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
plot(plot)
ggsave("figure_S2D.png", plot=plot, width = 6, height = 6, units="in") # CHANGE
dev.off()

