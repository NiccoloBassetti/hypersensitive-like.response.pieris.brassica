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




#### Statistical analysis ####
#____________________________#
packages=c("multcomp", "multcompView", "FSA", "rcompanion")
lapply(packages, library, character.only = TRUE)

mydata <- read.table("figure_3_dataset_PR1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_ICS1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_PR2.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_MYC2.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_VSP1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_VSP2.txt",header=TRUE,sep="")

tail(mydata)
names(mydata)

dim(mydata)
str(mydata)


# Filter data
data = as_tibble(mydata)

data = filter(data, Treatment=="Eggwash")
str(data)

# alternative model - separate time points
data_temp  =  as_tibble(mydata)
data_temp = filter(data_temp, Treatment=="Eggwash")

# Time point - 6 hours
data = filter(data_temp, Time_point=="6")
# Time point - 24 hours
data = filter(data_temp, Time_point=="24")
# Time point - 48 hours
data = filter(data_temp, Time_point=="48")


# Transform variables - TO DO after filtering dataset
data$Genotype = as.factor(data$Genotype)
data$Time_point = as.factor(data$Time_point)
data$Relative_exp = as.numeric(data$Relative_exp)
str(data)


# Check normality of the mydata
qqnorm(log2(data$Relative_exp))
qqline(log2(data$Relative_exp))
shapiro.test(log2(data$Relative_exp))		# test for normality (small P = NOT NORMAL)
fligner.test(log2(data$Relative_exp)~data$Genotype)	# test for equal variances (small P = NOT EQUAL VARIANCE)


# linear model    
lm1 <- lm(log2(Relative_exp) ~ Genotype * Time_point, data)
summary(lm1)
anova(lm1)
AIC(lm1) # goodness-of-fit measure - smaller values are better
BIC(lm1) # goodness-of-fit measure - smaller values are better

lm1 <- lm(log2(Relative_exp) ~ Genotype + Time_point, data)
summary(lm1)
anova(lm1)
AIC(lm1) # goodness-of-fit measure - smaller values are better
BIC(lm1) # goodness-of-fit measure - smaller values are better\

lm1 <- lm(log2(Relative_exp) ~ Genotype, data)
summary(lm1)
anova(lm1)
AIC(lm1) # goodness-of-fit measure - smaller values are better
BIC(lm1) # goodness-of-fit measure - smaller values are better\


coef(lm1) # coefficients of the model
confint(lm1) # confidence interval

# multiple comparison - 1 factor
pairs <- glht(lm1, linfct = mcp(Genotype = "Tukey"))
cld(pairs, level=0.05)
summary(pairs) 
confint(pairs)

# multiple comparison - 2 factors
library(lsmeans)
lsmeans = lsmeans(lm1, pairwise ~ Genotype + Time_point) # "+" means only additive effect and interaction effect
lsmeans  

lsmeans = lsmeans(lm1, pairwise ~ Genotype | Time_point) # "|" means only additive effect
lsmeans  

lsmeans = lsmeans(lm1, pairwise ~ Genotype)
lsmeans

# multiple comparison - 2 factors
cld(lsmeans, 
    alpa = 0.05,
    letters=letters,
    adjust  = "tukey")

# t-tests for multiple comparisons 
mydata <- read.table("figure_3_dataset_PR1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_ICS1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_PR2.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_MYC2.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_VSP1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_VSP2.txt",header=TRUE,sep="")

data  <-  as_tibble(mydata)

data = filter(data, Treatment=="Eggwash")


# Filter dataset to select timepoint to be tested
data = filter(data, Time_point=="6")
data = filter(data, Time_point=="24")
data = filter(data, Time_point=="48")

# Transform variables
data$Genotype = as.factor(data$Genotype)
data$Time_point = as.factor(data$Time_point)
data$Relative_exp = as.numeric(data$Relative_exp)
str(data)


t.test(log2(data$Relative_exp) ~ data$Genotype)



#### Plot results ####
#____________________#

packages = c("ggplot2", "svglite","agricolae","ggtext")
lapply(packages, library, character.only = TRUE)


mydata <- read.table("figure_3_dataset_PR1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_ICS1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_PR2.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_MYC2.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_VSP1.txt",header=TRUE,sep="")
mydata <- read.table("figure_3_dataset_VSP2.txt",header=TRUE,sep="")


# Reformat dataframe for plotting purpose
data_temp  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
data_temp$Relative_exp = as.numeric(data_temp$Relative_exp) # helps to plot x-axis variable in correct order
data_temp$Genotype = factor(data_temp$Genotype, levels = unique(data_temp$Genotype)) # helps to plot x-axis variable in correct order
data_temp$Time_point = factor(data_temp$Time_point, levels = unique(data_temp$Time_point)) # helps to plot x-axis variable in correct order
data_temp$Treatment = factor(data_temp$Treatment, levels = unique(data_temp$Treatment)) # helps to plot x-axis variable in correct order
data_temp
# data <- data %>% filter(`staining_type`=="TB")
# data <- data %>% dplyr::select(Species, staining, HR, staining) %>% gather(HR, staining, key = class, value = value)

# Filter data to keep only Eggs or Egg wash treatment
head(data_temp)
tail(data_temp)

trait = "Eggwash"

data = filter(data_temp, Treatment == trait)


head(data, 20)
tail(data)

str(data)

# sort dataset
#data = data %>% arrange(desc(treatment))


### Boxplot

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
  geom_boxplot(aes(x=Time_point, y=Relative_exp
                         #,fill=Genotype # ???
                         ,color=Genotype
        ),
        #colour="black", # change to black if you need to plot >1 Treatments
        outlier.colour=NA) +
  geom_point(mapping=aes(x=Time_point,y=Relative_exp
                         ,fill=Genotype
                         #,color=Genotype
                          ),
              colour="black", # black outline for points
              size=2,
              shape=19,
              alpha=0.5,
              position=position_jitterdodge(jitter.width=0.15,jitter.height=0,dodge.width=0.75,seed=NA)) +
  labs(x="", y="Relative  expression", size=size) +
  #scale_fill_manual(values=rep("white",length(levels(data$Time_point)))) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
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

ggsave("PLOT_PR1_FEEsubmission.tiff", plot=plot, width = 6, height = 6, units="in") # CHANGE
ggsave("PLOT_ICS1_FEEsubmission.tiff", plot=plot, width = 6, height = 6, units="in") # CHANGE
ggsave("PLOT_PR2_FEEsubmission.tiff", plot=plot, width = 6, height = 6, units="in") # CHANGE
ggsave("PLOT_MYC2_FEEsubmission.tiff", plot=plot, width = 6, height = 6, units="in") # CHANGE
ggsave("PLOT_VSP1_FEEsubmission.tiff", plot=plot, width = 6, height = 6, units="in") # CHANGE
ggsave("PLOT_VSP2_FEEsubmission.tiff", plot=plot, width = 6, height = 6, units="in") # CHANGE

dev.off()
