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
  temp$genotype <- factor((strsplit(file_list[i], split = "_")[[1]])[1])
  
  colnames(temp)=c("param", "plant", "value", "genotype")
  data_table = rbind(data_table, temp)
}

rm(temp) #remove data frame created in for loop

data_npq <- subset(data_table, data_table$param=="NPQ_Lss")[,3:4]

#Linear Discriminant Analysis for classification
#divide data frame into train and test groups
train_idx <- sample(nrow(data_npq), 0.8*nrow(data_npq))
train <- data_npq[train_idx,]
test <- data_npq[-train_idx,]

#build model (note the probabilities are equal for each group)
lda.fit <- lda(genotype~value, data = data_npq, subset = train_idx)

#prediction
lda.pred <- predict(lda.fit, test)
lda.class <- lda.pred$class
table(lda.class, test$genotype) #confusion matrix
