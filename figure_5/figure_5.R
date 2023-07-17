############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "RawdataFigure5.csv"
###
### Author R script:   Dr. Lotte Caarls - lotte.caarls@wur.nl 
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


# import dataset
file <- "RawdataFigure5.csv"
data <- read.csv(file,sep=";")

#rename 
data <- data  %>%  
  rename(experiment=ï..experiment)
#give me a summary of the data
summary(data)
data$treatment <- as.factor(data$treatment)
data$experiment<- as.factor(da
                            ta$experiment)
#make an occurance (HR yes) column
data$HR <- as.factor(ifelse(data$score>=2,"1","0"))
data$noHR <- as.factor(ifelse(data$score<2,"1","0"))


##Figure 5A - wash eggs of different age and hatched eggs###
#select experiment for figure 5A
data1 <- subset(data, experiment=="days")
levels(data1$treatment)

#remove other levels as factor
data1$treatment <- factor(data1$treatment,levels = c("1 day old","2 days old","3 days old","4 days old","5 days old","hatched"))

#summarize information for species and store in 'supdata': N, N of plants, HRy, HR ratio, mean, SD and SE
supdata <- data1 %>%
  group_by(treatment) %>%
  summarize(N    = length(score),
            HRy   = length(HR[HR == 1]),
            HRn   = length(HR[HR == 0]),
            HRratio = HRy/N,
            Mean = mean(score),
            SD = sd(score),
            SE = SD/N)
as.data.frame(supdata)
#### Statistical analysis ####
packages=c("multcomp", "multcompView", "FSA", "rcompanion")

lapply(packages, library, character.only = TRUE)
# Check normality of the data
qqnorm(data1$score)
qqline(data1$score)
shapiro.test(data1$score)		# test for normality (small P = NOT NORMAL)
fligner.test(data1$score~data1$treatment)	# test for equal variances (small P = NOT EQUAL VARIANCE)

#non-parametric test
kruskal.test(data1$score~data1$treatment)

# pairwise comparison: Withney-Mann U test/Wilcoxon rank sum test
PT = pairwise.wilcox.test(data1$score, data1$treatment,
                          p.adjust.method = "BH")
PT = PT$p.value
PT1 = round(fullPTable(PT),3) # require "rcompanion"

multcompLetters(PT1,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE) # require "multcompView"
PT1

#produce graph
packagesgraph=c("ggplot2", "reshape2", "RColorBrewer")

lapply(packagesgraph, library, character.only = TRUE)

#first prepare table
#decide which variables to use in the table
rowvariable <- data1$score
colvariable <- data1$treatment

# First make a table to count the scores per variable
table1 <- table(rowvariable,colvariable)
table1
#make a proportions table,  the second argument is if totals should come from rows or columsn 
prop <- prop.table(table1,2)
perc <- prop*100

#Make  barplot with ggpplot
#Need to convert the table into a frame
data2 <-as.data.frame.matrix(perc)
#add score as seperate variable
data2$score<-rownames(data2)
#convert to long format
data2long <- melt(data2, id.vars=c("score"), value.name = "proportion")
#paste variable as name column 2
names(data2long)[2] <-paste("treatment")
# make the plot
# the order of the factors is the "wrong" way around per default, so reverse is TRUE
display.brewer.pal(n=9,"Blues")
my.blue <- brewer.pal(n=9,"Blues")[c(2,3,5,7,9)]
bp <- ggplot() + geom_bar(aes(y = proportion,
                              x = treatment,
                              fill = score),
                          data=data2long,
                          stat="identity",
                          position=position_stack(reverse=TRUE), width = 0.6) + # change the width of the bars
  scale_fill_manual(values=my.blue,name = 'class', guide = guide_legend(reverse=TRUE)) + # change the colours and reverse the legend
  scale_y_continuous(name = "% plants in class", expand = c(0,0), limits = c(0,101)) + 
  scale_x_discrete(name = "")
bp 

require(svglite)

# save plot as image, set measurements
graphname <- "Caarls1_days.svg"
ggsave(graphname, width = 300, height = 150, units = "mm")
####################################################
##Figure 5B - wash eggs, ARG, ovarian eggs, BC ###
#####################################################
data3 <- subset(data, experiment=="ARG dissection")
levels(data3$treatment)
#remove others as factors
data3$treatment <- factor(data3$treatment,levels = c("egg wash","ARG","eggs ovary","Bursa Cop","ARG unmated"))

# Check normality of the data
qqnorm(data3$score)
qqline(data3$score)
shapiro.test(data3$score)		# test for normality (small P = NOT NORMAL)
fligner.test(data3$score~data3$treatment)	# test for equal variances (small P = NOT EQUAL VARIANCE)

#non-parametric test
kruskal.test(data3$score~data3$treatment)


# pairwise comparison: Withney-Mann U test/Wilcoxon rank sum test
PT = pairwise.wilcox.test(data3$score, data3$treatment,
                          p.adjust.method = "BH")
PT = PT$p.value
PT1 = round(fullPTable(PT),3) # require "rcompanion"

multcompLetters(PT1,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE) # require "multcompView"
PT1

##prepare graph Figure 5B
#first prepare table
#decide which variables to use in the table
rowvariable <- data3$score
colvariable <- data3$treatment

# First make a table to count the scores per variable
table1 <- table(rowvariable,colvariable)
table1
#make a proportions table,  the second argument is if totals should come from rows or columsn 
prop <- prop.table(table1,2)
perc <- prop*100

#I need to convert the table into a frame
data2 <-as.data.frame.matrix(perc)
#add score as seperate variable
data2$score<-rownames(data2)
#convert to long format
data2long <- melt(data2, id.vars=c("score"), value.name = "proportion")

#paste variable as name column 2
names(data2long)[2] <-paste("treatment")
# make the plot
# the order of the factors is the "wrong" way around per default, so reverse is TRUE
Fig5B <- ggplot() + geom_bar(aes(y = proportion,
                              x = treatment,
                              fill = score),
                          data=data2long,
                          stat="identity",
                          position=position_stack(reverse=TRUE), width = 0.6) + # change the width of the bars
  scale_fill_manual(values=my.blue,name = 'class', guide = guide_legend(reverse=TRUE)) + # change the colours and reverse the legend
  scale_y_continuous(name = "% plants in class", expand = c(0,0), limits = c(0,101)) + 
  scale_x_discrete(name = "")
Fig5B 
# save plot as image, set measurements
graphname <- "Caarls1_ARG.svg"
ggsave(graphname, width = 300, height = 150, units = "mm")

##################################################
##Figure 5C - wash of mated vs virgin eggs  ###
#####################################################
data4 <- subset(data, experiment=="virgin egg wash")
levels(data4$treatment)
#remove those as factors
data4$treatment <- factor(data4$treatment,levels = c("egg wash","eggs unfertilized"))

# Check normality of the data
qqnorm(data4$score)
qqline(data4$score)
shapiro.test(data$score)		# test for normality (small P = NOT NORMAL)
fligner.test(data4$score~data4$treatment)	# test for equal variances (small P = NOT EQUAL VARIANCE)

#non-parametric test
kruskal.test(data4$score~data4$treatment)


# pairwise comparison: Withney-Mann U test/Wilcoxon rank sum test
PT = pairwise.wilcox.test(data4$score, data4$treatment,
                          p.adjust.method = "BH")
PT = PT$p.value
PT1 = round(fullPTable(PT),3) # require "rcompanion"

multcompLetters(PT1,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE) # require "multcompView"
PT1

#decide which variables to use in the table
rowvariable <- data4$score
colvariable <- data4$treatment
# First make a table to count the scores per variable
table1 <- table(rowvariable,colvariable)

#make a proportions table,  the second argument is if totals should come from rows or columsn 
prop <- prop.table(table1,2)
perc <- prop*100

#I need to convert the table into a frame
data2 <-as.data.frame.matrix(perc)
#add score as seperate variable
data2$score<-rownames(data2)
#convert to long format
data2long <- melt(data2, id.vars=c("score"), value.name = "proportion")

#paste variable as name column 2
names(data2long)[2] <-paste("treatment")
# make the plot
# the order of the factors is the "wrong" way around per default, so reverse is TRUE
Fig5C <- ggplot() + geom_bar(aes(y = proportion,
                                 x = treatment,
                                 fill = score),
                             data=data2long,
                             stat="identity",
                             position=position_stack(reverse=TRUE), width = 0.6) + # change the width of the bars
  scale_fill_manual(values=my.blue,name = 'class', guide = guide_legend(reverse=TRUE)) + # change the colours and reverse the legend
  scale_y_continuous(name = "% plants in class", expand = c(0,0), limits = c(0,101)) + 
  scale_x_discrete(name = "")
Fig5C 
# save plot as image, set measurements
graphname <- "Fig5C.svg"
ggsave(graphname, width = 300, height = 150, units = "mm")

############################################
##Figure 5D - wash glue and glue removed  ###
#####################################################

data5 <- subset(data, experiment=="glue")
levels(data5$treatment)
#remove others as factors
#change order
data5$treatment <- factor(data5$treatment,levels = c("egg wash","wash glue","glue removed (bleach)","glue removed (NaPO4)"))
data5
# Check normality of the data
qqnorm(data5$score)
qqline(data5$score)
shapiro.test(data$score)		# test for normality (small P = NOT NORMAL)
fligner.test(data5$score~data5$treatment)	# test for equal variances (small P = NOT EQUAL VARIANCE)

#non-parametric test
kruskal.test(data5$score~data5$treatment)


# pairwise comparison: Withney-Mann U test/Wilcoxon rank sum test
PT = pairwise.wilcox.test(data5$score, data5$treatment,
                          p.adjust.method = "BH")
PT = PT$p.value
PT1 = round(fullPTable(PT),3) # require "rcompanion"

multcompLetters(PT1,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE) # require "multcompView"
PT1
#decide which variables to use in the table
rowvariable <- data5$score
colvariable <- data5$treatment
# First make a table to count the scores per variable
table1 <- table(rowvariable,colvariable)

#make a proportions table,  the second argument is if totals should come from rows or columsn 
prop <- prop.table(table1,2)
perc <- prop*100

#I need to convert the table into a frame
data2 <-as.data.frame.matrix(perc)
#add score as seperate variable
data2$score<-rownames(data2)
#convert to long format
data2long <- melt(data2, id.vars=c("score"), value.name = "proportion")

#paste variable as name column 2
names(data2long)[2] <-paste("treatment")
# make the plot
# the order of the factors is the "wrong" way around per default, so reverse is TRUE
Fig5D <- ggplot() + geom_bar(aes(y = proportion,
                                 x = treatment,
                                 fill = score),
                             data=data2long,
                             stat="identity",
                             position=position_stack(reverse=TRUE), width = 0.6) + # change the width of the bars
  scale_fill_manual(values=my.blue,name = 'class', guide = guide_legend(reverse=TRUE)) + # change the colours and reverse the legend
  scale_y_continuous(name = "% plants in class", expand = c(0,0), limits = c(0,101)) + 
  scale_x_discrete(name = "")
Fig5D 
# save plot as image, set measurements
graphname <- "Fig5D.svg"
ggsave(graphname, width = 300, height = 150, units = "mm")
#################################################################

