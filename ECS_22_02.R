setwd("C:/Users/maciek/Desktop/Studia_dok/DualPAM/ECS/experiment_npq")   #source files location

library(ggplot2)
library(RColorBrewer)
library(stringr)
library(grid)
library(plyr)
library(reshape2)
library(viridis)

file_list <- list.files(pattern=".CSV")   #load data files

data_f = data.frame(   #create data frame to store data 
  "Time" = numeric(),   
  "ECS" = numeric(),   #relative absorbance values
  "genotype" = factor(),   #plant genotype
  "al"= factor(),   #actinic light intensity during measurement
  "trt"= factor(),   #Control, Excess light (BL, RL), days after UV trt
  "read" = factor()   #plant number
)
for(i in 1:length(file_list)){
  temp = as.data.frame(read.table(file_list[i], header=F, sep=";", skip=2))[,c(1,2)]
  
  time <- temp$V1[seq(1,length(temp$V1),10)]   #every value as a mean from 10 time points to reduce noise
  #time <- time[-1]
  matr <- matrix(temp$V2,byrow=T,ncol=10)   
  #matr <- matr[-nrow(matr),]
  data_temp <- rowMeans(matr)
  data_temp <- data.frame(cbind(time,data_temp))
  
  data_temp$V3 = (strsplit(file_list[i], split =c("_|\\."))[[1]])[1]
  data_temp$V4 = (strsplit(file_list[i], split =c("_|\\."))[[1]])[3]
  data_temp$V5 = (strsplit(file_list[i], split =c("_|\\."))[[1]])[2]
  data_temp$V6 = (strsplit(file_list[i], split =c("_|\\."))[[1]])[4]
  colnames(data_temp) = c("Time","ECS", "genotype", "al", "trt", "read")
  
  data_f = rbind(data_f,data_temp)   #bind data to main data frame
}

rm(matr,data_temp,temp,time)   #remove temporary data frames

data_f$genotype <- factor(data_f$genotype, levels = c('col0','npq4','oePsbS'))
data_f$al <- factor(data_f$al, levels = c("AL200", "AL600"))
data_f$trt <- factor(data_f$trt)
data_f$read <- factor(data_f$read)

#Further data manipulation

data_ECS_kinetics <- subset(data_f, data_f$Time > 130 & data_f$Time < 165) #subset only data interesting for analysis

baseline_start <- subset(data_ECS_kinetics, data_ECS_kinetics$Time < 139)   #values before switching light off
baseline_start <- ddply(baseline_start,~genotype+al+trt+read,summarise,mean=mean(ECS))   #calculate means of baselines (unsteady ECS curve)

relax_min <- ddply(data_ECS_kinetics[data_ECS_kinetics$Time < 142,],
                   ~genotype+al+trt+read,summarise, ECS=min(ECS))   #curve decay after switching light off -- minimal value

baseline_dark <- subset(data_ECS_kinetics, data_ECS_kinetics$Time > 155)
baseline_dark <- ddply(baseline_dark, ~genotype+al+trt+read, summarise, ECS=mean(ECS))   #steady value of dark relaxed ECS

pmf <- baseline_start[,1:3]   #create data frame with proton motive force components
pmf$whole <- baseline_start$mean - relax_min$ECS   #whole pmf (sum of potential and pH, from baseline in light to minimal value during dark relaxation)
pmf$pH <- baseline_dark$ECS - relax_min$ECS
pmf$pH_rate <- pmf$pH/pmf$whole   #ratio of pH part of pmf to the whole value of pmf (to allow comparison between measurements)


#Ggplot theme

mytheme <- theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color="black"),
    plot.margin = unit(c(0,0,0,0),"mm"),
    #legend.title = element_blank(),
    #legend.position = c(0.9,0.9),
    strip.text.y = element_blank(),
    #strip.text.x = element_blank(),
    strip.background = element_blank()
    #legend.position = "none"
  )

#Plots
pd = position_dodge(width = .9)
gg_ECS <- ggplot(pmf, aes(x=trt, y=pH_rate, fill=genotype)) + geom_boxplot(position = pd) 
gg_ECS + facet_wrap(~al) + mytheme + scale_fill_viridis(discrete = TRUE)


