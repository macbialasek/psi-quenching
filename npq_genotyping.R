setwd("C:/Users/maciek/Desktop/Studia_dok/Barczak/quenching/Second_iteration/AL600")
library(MASS) 
library(reshape2)
library(randomForest)
library(caret)
library(e1071)

file_list <- list.files(pattern = c("Control"))  #read files; control measuremnets

data_table = data.frame(   #create data frame for data storage
  "param" = factor(),
  "plant"= factor(),
  "value"= numeric(),
  "genotype"=factor(),
  check.names= FALSE
)

for(i in 1:length(file_list)){  #bind data from files to data frame
  
  temp <- as.data.frame(read.table(file_list[i], skip=3, header=F))   #temp for temporary data storage
  temp <- rbind(temp[42,],temp[47,],temp[52,], temp[64,])   #read only four parameters (npq, psii, fv/fm, fv'/fm')
  temp <- melt(temp, id.vars=("V1"))   #wide to long type of table
  temp$genotype <- factor((strsplit(file_list[i], split = "_")[[1]])[1])   #names of genotypes from file name
  
  colnames(temp)=c("param", "plant", "value", "genotype")
  data_table = rbind(data_table, temp)   #bind temporary data frame with final df
}

rm(temp) #remove data frame created in for loop

data_npq <- subset(data_table, data_table$param=="NPQ_Lss")[,3:4]   #subset data frame; only npq values

#Linear Discriminant Analysis for classification
#data spliting method from caret
set.seed(3456)   #for reproducibility
trainIndex <- createDataPartition(data_npq$genotype, p = .8,   #split df to train and test subset
                                  list = FALSE,
                                  times = 1)
trainCV <- createMultiFolds(data_npq$genotype)   #for cross-validation

#build model (note the probabilities are equal for each group)
lda.fit <- lda(genotype~value, data = data_npq, subset = trainIndex)

lda.pred <- predict(lda.fit, data_npq[-trainIndex,])   #test model
lda.class <- lda.pred$class   #I need class prediction, not probabilities
confusionMatrix(lda.class, data_npq[-trainIndex,]$genotype) #confusion matrix

#LDA model with caret syntax
ctrl <- trainControl(method = "cv",
                     classProbs = TRUE,
                     summaryFunction = multiClassSummary)

ldaModel <- train(genotype ~ value, data = data_npq[trainIndex,],
                  method = "lda",
                  metric = "ROC",
                  trControl = ctrl)

getTrainPerf(ldaModel)
ldaPred <- predict(ldaModel,data_npq[-trainIndex,])
confusionMatrix(ldaPred, data_npq[-trainIndex,]$genotype)

#random forest method
set.seed(1)
rf.npq <- randomForest(genotype ~ value, data = data_npq, subset = trainIndex, 
                       importance = TRUE, ntree = 1000)
yhat.rf <- predict(rf.npq, newdata = data_npq[-trainIndex,])
confusionMatrix(yhat.rf,data_npq[-trainIndex,]$genotype)

#random forest caret syntax
ctrl <- trainControl(method = "cv",
                     classProbs = TRUE,
                     summaryFunction = multiClassSummary)

rfModel <- train(genotype ~ value, data = data_npq[trainIndex,],
                  method = "rf",
                  metric = "ROC",
                  trControl = ctrl)

getTrainPerf(rfModel)
rfPred <- predict(rfModel,data_npq[-trainIndex,])
confusionMatrix(rfPred, data_npq[-trainIndex,]$genotype)