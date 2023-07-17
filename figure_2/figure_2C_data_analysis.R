############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "figure_2C_dataset.txt"
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

mydata <- read.table(file = "figure_2C_dataset.txt",sep = "\t", header = TRUE) # !!! test sep character

str(mydata)
head(mydata)


# Reformat dataframe for plotting purpose
data_temp  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
data_temp$Timepoint = factor(data_temp$Timepoint, levels = unique(data_temp$Timepoint)) # helps to plot x-axis variable in correct order
data_temp$Treatment = factor(data_temp$Treatment, levels = unique(data_temp$Treatment)) # helps to plot x-axis variable in correct order
data_temp


# Reformat dataset - insert a copy of Controls and label each copy with either as treated wiht "Eggs" or "Egg wash"
# Controls were of course untreated plants. The treatment labelling is needed to perform statistical testing. 
head(data_temp)
tail(data_temp)


data_eggs = filter(data_temp, Treatment == "eggs")
data_egg_wash = filter(data_temp, Treatment == "egg_wash")
data = rbind(data_eggs, data_egg_wash)
data = rbind(filter(filter(data_temp, Treatment == "Control"),Timepoint == "0"),data)  
data = rbind(filter(filter(data_temp, Treatment == "Control"),Timepoint == "0"),data)
data[1:4,5] = "eggs" 
data[5:8,5] = "egg_wash" 

head(data, 20)
tail(data)

str(data)


#### Statistical analysis ####
#____________________________#

packages=c("multcomp", "multcompView", "FSA", "rcompanion")
lapply(packages, library, character.only = TRUE)


# Check normality of the mydata
qqnorm(log2(data$Relative_expression))
qqline(log2(data$Relative_expression))
shapiro.test(log2(data$Relative_expression))		# test for normality (small P = NOT NORMAL)
fligner.test(log2(data$Relative_expression)~data$Timepoint)	# test for equal variances (small P = NOT EQUAL VARIANCE)


# Transform class variables
sapply(FUN = class, data)
data$Timepoint = factor(data$Timepoint, levels = unique(data$Timepoint))
str(data)


# linear model
lm1 <- lm(log2(Relative_expression) ~ Timepoint*Treatment, data)

summary(lm1)
anova(lm1)
AIC(lm1) # goodness-of-fit measure - smaller values are better
BIC(lm1) # goodness-of-fit measure - smaller values are better
coef(lm1) # coefficients of the model
confint(lm1) # confidence interval

# multiple comparison
pairs <- glht(lm1, linfct = mcp(Timepoint = "Tukey"))
cld(pairs, level=0.05)
summary(pairs) 
confint(pairs)

# CONCLUSION
# Data are normally distribute and can be analysed with two-way ANOVA.
# Only Timepoint have a significant effect. Treatments can then analysed separately.
# A comparison between treatments is also possible but it was not investigated. 


#### Statistical analysis - separate treatments ####
#____________________________#

# Filter data to keep only Eggs or Egg wash treatment
head(data)
tail(data)


data_2 = filter(data, Treatment == "eggs")

  data_2 = filter(data, Treatment == "egg_wash") # repeat analysis selecting "egg wash"

head(data_2, 20)
tail(data_2)

str(data_2)

# Test for normality and homogeneity of variances
qqnorm(log2(data_2$Relative_expression))
qqline(log2(data_2$Relative_expression))
shapiro.test(log2(data_2$Relative_expression))		# test for normality (small P = NOT NORMAL)
fligner.test(log2(data_2$Relative_expression)~data_2$Timepoint)	# test for equal variances (small P = NOT EQUAL VARIANCE)


# CONCLUSION: Data looks normal after transformation, nevertheless were analysed with non-parametric tests
# because the few data points for every treatment


# non-parametric test
kruskal.test(data_2$Relative_expression~data_2$Timepoint)

# pairwise comparison 1: Dunn test
PT = dunnTest(data_2$Relative_expression~data_2$Timepoint,
              data=data_2,
              method="bh")    # require "FSA", p.adjust to be used to adjust p-values
PT
PT = PT$res

# CONCLUSION
# EGG: for each timepoint, only comparison against control (0h) were considered relevant for the study
# EGG WASH: for each timepoint, only comparison against control (0h) were considered relevant for the study



#### Plot results ####
#____________________#

packages = c("ggplot2", "svglite","agricolae","ggtext")
lapply(packages, library, character.only = TRUE)

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
                   #,fill=Treatment # ???
                   ,color=Treatment
                   ),
               #colour="black", # change to black if you need to plot >1 Treatments
               outlier.colour=NA) +
  geom_point(mapping=aes(x=Timepoint,y=Relative_expression
                         ,fill=Treatment
                         #,color=Genotype
                         ),
             colour="black", # black outline for points
             size=2,
             shape=19,
             alpha=0.5,
             position=position_jitterdodge(jitter.width=0.15,jitter.height=0,dodge.width=0.75,seed=NA)) +
  labs(x="", y="Relative  expression", size=size) +
  #scale_fill_manual(values=rep("white",length(levels(data$Timepoint)))) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_fill_brewer(palette=cbPalette) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_colour_brewer(palette=cbPalette) + # set color for line/points - 
  scale_colour_manual(values = c("#E31A1C", "#1F78B4")) + # set color for line/points -
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
  #theme(legend.position="none") + # change to set position, "none" remove legend
  theme(legend.position = c(0.30, 0.90)) + # remove this to override "none"
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size=(size), face="plain")) + # change "face" value for style text
  #scale_x_discrete(limits=c("0h","3h","6h","24h","48h")) + # reorder scale x and legend
  #scale_y_continuous(#breaks=c(0,seq(100,300, by=100)),
  #breaks=c(seq(0,200, by=50)),
  #limits=c(0,300)
  #)+
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
plot(plot)
ggsave("figure_2C.png", plot=plot, width = 6, height = 6, units="in")
dev.off()

