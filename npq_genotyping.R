setwd("C:/Users/maciek/Desktop/Studia_dok/Barczak/quenching/Second_iteration/AL600")
library(MASS)
library(reshape2)

#read files with results from control measurements

file_list <- list.files(pattern = c("Control"))

data_table = data.frame(
  "param" = factor(),
  "plant"= factor(),
  "value"= numeric(),
  "genotype"=factor(),
  check.names= FALSE
)

for(i in 1:length(file_list)){
  temp <- as.data.frame(read.table(file_list[i], skip=3, header=F))
  temp <- rbind(temp[42,],temp[47,],temp[52,], temp[64,])
  temp <- melt(temp, id.vars=("V1"))
  temp$genotype <- (strsplit(file_list[i], split = "_")[[1]])[1]
  
  colnames(temp)=c("param", "plant", "value", "genotype")
  data_table = rbind(data_table, temp)
}

rm(temp)
