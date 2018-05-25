# Setting up directory
getwd()
setwd("D:/DS/insa/Os meus documentos/Google Drive/R")
setwd("C:/Users/Asif/Desktop/R")
setwd("C:/Users/Asif/Google Drive/R")
# load libraries
library(caret) # contain many model desiging functions
library(pROC)  # function to analyze matrix to calcualte the AUC(area under the curve?) 
library(rattle) # to draw nice tree
library(AppliedPredictiveModeling)
library(psych) 
library(PerformanceAnalytics)
library(dplyr)
library(mlr)
library(unbalanced)
library(rpart.plot)
library(rpart)
library(ROCR)
#library(installr)
#updateR()
#version
# read data
AGP_data<-read.table("test.txt", header = T, sep = "\t")

# to know about data
dim(AGP_data)
str(AGP_data)
sapply(AGP_data, class)
summary(AGP_data)
names(AGP_data)

# class attribute statistics
table(AGP_data[ncol(AGP_data)])
prop.table(table(AGP_data[ncol(AGP_data)]))

######### make the family type a numeric variable
table(AGP_data$fam_type)
AGP_data$fam_type <- ifelse(AGP_data$fam_type=='MPX',1,ifelse(AGP_data$fam_type=='SPX',2, 
                                                              ifelse(AGP_data$fam_type=='UKN',3,4)))

names(AGP_data)
######################

##############################################
# Model Designing
###############################################
###########################################

# get names of all caret supported models 
names(getModelInfo())

#titanicDF$Survived <- ifelse(titanicDF$Survived==1,'yes','nope')
#table(titanicDF$Survived)
# pick model gbm and find out what type of model it is
getModelInfo()$gbm$type

AGP_data<-read.table("test.txt", header = T, sep = "\t")
# to know about data
dim(AGP_data)
str(AGP_data)
sapply(AGP_data, class)
summary(AGP_data)
names(AGP_data)

# class attribute statistics
table(AGP_data[ncol(AGP_data)])
prop.table(table(AGP_data[ncol(AGP_data)]))


# correlation among predicters
#Drawing the correlation among all variables
CorPlot<-pairs.panels(AGP_data[1:8]) # select columns 1-4
##########
png(filename="CorPLot.png")
pairs.panels(AGP_data[1:8])
dev.off()


################


outcomeName<-names(AGP_data[ncol(AGP_data)])
predictorsNames <- names(AGP_data)[names(AGP_data) != outcomeName]
str(predictorsNames)

AGP_data<-read.table("test.txt", header = T, sep = "\t")
#rounding the data frame
AGP_data<-data.frame(lapply(AGP_data, function(y) if(is.numeric(y)) round(y, 1) else y)) 

#get 10% data for testing
si<-round(.1*nrow(AGP_data),0)
###############################
asd <- sample(seq_len(nrow(AGP_data)), size = si, replace = F, prob=AGP_data$NervousSystem) 


testing<-AGP_data[asd,]
str(testing)
train<-AGP_data[-asd,]
str(train)
train$NervousSystem<- ifelse(train$NervousSystem=="yes",1 , 0)
#I have to reduce the positive class
training<-OverSam(train)

rownames(training)<-NULL
training$NervousSystem<- ifelse(training$NervousSystem==1,"yes","no")

str(training)
table(training$NervousSystem)
inc<-0.5
x = list(1, 2, 3, 3,2)
aucValues<-c()
Results<-data.frame ()
percentData<-c()

for(i in 1:5){
  sii<-round(inc*nrow(training),0)
  #I can add probabilty arugument
  IDs <- sample(seq_len(nrow(training)), size = sii, replace = F) 
  
  
  trainingSub<-training[IDs,]
  
  model<-rpart(NervousSystem~., data=trainingSub, method="class")
  

  #Evaluation
  # Calculating confusion matrix
  EvaTest<-predict(model, testing,type="class")
  table(testing$NervousSystem, EvaTest)
  xtab<-table(testing$NervousSystem, EvaTest)
  #Function call
  RPartlist<-ConMat(xtab)
  Results <- rbind(Results, (unlist(RPartlist)))
  #       names(testing)
  ## calculating ROC curves
  #plot(rp)
  q<-predict(model,  newdata=testing,type = "prob")
  pred <- prediction(q[,2], testing$NervousSystem)
  str(pred)
  x[[i]]<-performance(pred, "tpr", "fpr")
  aucc<-performance(pred, "auc")
  aucValues<-rbind(aucValues,aucc@y.values[[1]])
  
  
  percentData<-rbind(percentData,(inc*100))
  
  
  
  inc<-inc+0.1
  
  
}
####################################
mean(aucValues)
colnames(aucValues)<-c("AUC")
colnames(Results) <- c("Accuracy", "Precision", "Recall", "Specificity", "F-measure")
colnames(percentData)<-c("DataPercentage")
Results2<-cbind(Results,aucValues,percentData)


Results3<-round(rbind(Results2,apply(Results2,2,mean)),3)
rownames(Results3) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5", "mean")

write.csv(Results3, "res.csv")
#

########################################
#plot results
names(Results3)
Results3<-Results3[nrow(Results3)-1,]
plot(Results3$DataPercentage,Results3$AUC, type="l")


##############################################################
# All functions
#####################################################################################################
########################################################################################################
#Random OverSampling
OverSam<-function(train){
  
  n<-ncol(train)
  output<-train[ ,n]
  input<-train[ ,-n]
  data<-ubOver(X=input, Y= output, k=0)
  str(data)
  newData<-cbind(data$X, data$Y)
  colnames(newData)<-c(names(train))
  names(newData)
  dim(newData)
  #  newData$NervousSystem<- ifelse(newData$NervousSystem==1,"yes","no")
  
  return(newData)
  
}

# Function to calculate confusion matrix
ConMat<-function(x){
  # Function created to calculate precision, recall and others stuff
  #Its uses confusionMatrix fucntion from CART package
  #take confusion table as input
  #returen a list which should be unlist before writting to data frame
  result <-confusionMatrix(x)
  attributes(result)
  precision <- result$byClass['Pos Pred Value']    
  recall <- result$byClass['Sensitivity']
  f_measure <- 2 * ((precision * recall) / (precision + recall))
  str(result)
  Spec<-result$byClass['Specificity']
  Accuarcy<- result$overall['Accuracy']
  ff<-list ( Accuarcy, precision, recall, Spec, f_measure)
  return (ff)
}
