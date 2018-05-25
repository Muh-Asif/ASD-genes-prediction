getwd()
setwd("D:/DS/insa/Os meus documentos/Google Drive/R")
setwd("C:/Users/Asif/Desktop/R")
setwd("C:/Users/Asif/Google Drive/R")
setwd("C:/Users/ASUS/Google Drive/R")
library(caret) # contain many model desiging functions
library(pROC)  # function to analyze matrix to calcualte the AUC(area under the curve?) 
library(rattle) # to draw nice tree
library(dplyr)
library(unbalanced)
library(rpart.plot)
library(rpart)
library(ROCR)
library(e1071)
library(klaR)
library(kernlab)
library(randomForest)
library(klaR)
library(RWeka)
library(psych)
library(tidyr)

library(rpart)
library(e1071)
library(tidyr)
library(randomForest)
library(caret)
library(mlr)
library(unbalanced)
#library(scales) 
library(ROCR)



dataIn<-read.table("88.txt", header = T, sep = "\t")
Genedata<-read.table("SimplexVsMulti.txt", header = T, sep = "\t")

Cate4and3<-Genedata[Genedata[,2]== 3 | Genedata[,2]== 4 | Genedata[,2]== "S",]

#Cate4and3<-Genedata[Genedata[,2]== 3 | Genedata[,2]== 4 ,]

#write.csv(Cate4and3,"159.CSV")
dim(Cate4and3)
names(Cate4and3)

Rem<-Cate4and3[,1]
length(Rem)
class(Rem)
dim(dataIn)
names(dataIn)
dataIn[1:6,1:6]

#res<-dataIn[ , !(names(dataIn) %in% paste("X",Rem,sep = ""))]
#dim(res)


#data<-res[!(res$ID  %in% Rem),]

# evaluating only on high confidence gened
data<-dataIn

dim(data)
#write.csv(data,"159.CSV")

data[1:6,1:6]

names(data)[2]<-paste("classVar")
########################################
names(data)
dim(data)
table(data[2])
table(data$classVar)
data[1:6,1:6]

# remove columns with all NAz
ind2 <- apply(data, 2, function(x) all(is.na(x)))
SemMatrixWithOutNa <- data[ ,!ind2 ]
dim(SemMatrixWithOutNa)
dat_full<-SemMatrixWithOutNa
# remove rows with any NAz
row.has.na <- apply(dat_full, 1, function(x){any(is.na(x))})
sum(row.has.na)
final.filtered <- dat_full[!row.has.na,]
dat_full<-final.filtered 
dim(dat_full)
#write.csv(dat_full, "159.CSV")

table(dat_full$classVar)

### feature selection
NumFold<-5
seedNum=341596324
#13120=0.649; 131205=0.659; 131239= 0.646; 17000000=0.662
#19000000=0.648; 3596841=0.659; 341596324=0.655
#?ubOver

#OverSampleRes<-RF_overSampleModel(NumFold,dat_full, seedNum)
#OverSampleRes

res<-Navie.bayes_model(NumFold,dat_full, seedNum)
res
rownames(res)<-c("Navie.Bayes.fold.1","fold.2","fold.3","fold.4","fold.5","Mean")

Rpart_res<-rpart_model(NumFold,dat_full, seedNum)
rownames(Rpart_res)<-c("CART.fold.1","fold.2","fold.3","fold.4","fold.5","Mean")
Rpart_res

#341596324: 0.727; 17000000=0.726, 131239=0.713; 90152584=0.727
#par(mfrow=c(2,2))
RF_res<-RF_model(NumFold,dat_full, seedNum)
RF_res
rownames(RF_res)<-c("RandomForest.fold.1","fold.2","fold.3","fold.4","fold.5","Mean")
RF_res

Radial_svm_res<-svm_model2(NumFold,dat_full, seedNum)
rownames(Radial_svm_res)<-c("Radial-SVM.fold.1","fold.2","fold.3","fold.4","fold.5","Mean")
Radial_svm_res

lin_svm_res<-linear_svm_model(NumFold,dat_full, seedNum)
rownames(lin_svm_res)<-c("Linear-SVM.fold.1","fold.2","fold.3","fold.4","fold.5","Mean")
lin_svm_res
write.csv(rbind(res, Rpart_res, RF_res,Radial_svm_res, lin_svm_res), "res.csv")


##################
############   functions
# SVM function
#########################
##################################
#RF_overSampleModel<-function(N_fold,AGP_data, SedVal){
  # AGP_data<-dat_full ####
  N_fold2<-N_fold
  outcomeName<-names(AGP_data[1])
  
  x = list() # to store ROC curve
  aucValues<-c() # to save auc value
  Results<-data.frame () # to save results from confusion
  set.seed(SedVal)#
  folds<-createFolds(AGP_data[,outcomeName], k = N_fold2, list = TRUE, returnTrain = FALSE)
  
  for (i in 1:N_fold2){
    
    # split data into training and testing chunks
    
    testing<- AGP_data[folds[[i]],]#i
    dim(testing)
    
    
    ################
    
    training<- AGP_data[-folds[[i]],]#i
    ####################### start of over sample
    train2<-training
    #call to overSampling method
    train2$classVar<- ifelse(train2$classVar=="ASD",1 , 0)
    #training<- UnderSam(train2)
    training<-OverSam(train2)
    
    rownames(training)<-NULL
    training$classVar<- ifelse(training$classVar==1,"ASD","nonASD")
    training$classVar<-as.factor(training$classVar)
    ###### end of over sample
    ###########
    #model desiging
    model<-randomForest(classVar~., data=training)
    summary(model)
    str(model)
    #fancyRpartPlot(model)
    #rpart.plot(model)
    
    #Evaluation
    # Calculating confusion matrix
    EvaTest<-predict(model, testing,type="class")#
    table(testing$classVar, EvaTest)
    xtab<-table(testing$classVar, EvaTest)
    #Function call
    RPartlist<-ConMat(xtab)
    Results <- rbind(Results, (unlist(RPartlist)))
    #       names(testing)
    ## calculating ROC curves
    #plot(rp)
    
    q<-predict(model,  newdata=testing,type = "prob")
    str(q)
    
    #    ROCRpred   <- ROCR::prediction(q[,2], testing$classVar)
    #   AUC        <- as.numeric(ROCR::performance(ROCRpred, "auc")@y.values)
    
    pred <- prediction(q[,2], testing$classVar)
    ##############
    
    ###############
    x[[i]]<-performance(pred, "tpr", "fpr")
    aucc<-performance(pred, "auc")
    aucValues<-rbind(aucValues,aucc@y.values[[1]])
    
    
    #ROC curve
    # png(filename="curve.png")
    
    if (i==1){
      plot(x[[i]], col=i)
    }
    else{
      plot(x[[i]], col=i,add=T)
      abline(0, 1, lty = 2)
      
    }
    #    dev.off()
    
    
  }
  
  
  colnames(aucValues)<-c("AUC")
  Results2<-cbind(Results,aucValues)
  Results3<-round(rbind(Results2,apply(Results2,2,mean)),3)
  colnames(Results3) <- c("Accuracy", "Precision", "Recall", "Specificity", "F-measure","AUC")#
  rownames(Results3) <- c(1:N_fold2, "mean")
  #  FileName<-paste(Pathway,"csv", sep=".")
  write.csv(Results3, file="res.CSV")
  return(Results3)
  
}# end of raprt function
#################

##################################
#OverSam<-function(train){
  #n<-ncol(train)
  output<-train[,1]
  input<-train[,-1]
  input[1:5,1:5]
  
  #  data<-ubOver(X=input, Y= output, k=0)
  data<-ubUnder(X=input, Y= output, perc = 30,  method = "percPos")
  str(data)
  
  newData<-cbind(data$Y,data$X)
  colnames(newData)<-c(names(train))
  #  newData$NervousSystem<- ifelse(newData$NervousSystem==1,"yes","no")
  
  return(newData)
  
  
}


#   rpart function
rpart_model<-function(N_fold,AGP_data, SedVal){
  # AGP_data<-NewData
  N_fold2<-N_fold
  outcomeName<-names(AGP_data[1])
  
  x = list() # to store ROC curve
  aucValues<-c() # to save auc value
  Results<-data.frame () # to save results from confusion
  set.seed(SedVal)
  folds<-createFolds(AGP_data[,outcomeName], k = N_fold2, list = TRUE, returnTrain = FALSE)
  
  for (i in 1:N_fold2){
    
    # split data into training and testing chunks
    
    testing<- AGP_data[folds[[i]],]
    dim(testing)
    testing<-testing[!(testing$ID  %in% Rem),]
    #write.csv(testing,"159.CSV")
    
    testing<-testing[,-1]
    ################
    
    training<- AGP_data[-folds[[i]],]
    training<-training[,-1]
    
    ###########
    #model desiging
    model<-rpart(classVar~., data=training, method="class")
    #fancyRpartPlot(model)
    #rpart.plot(model)
    
    
    #Evaluation
    # Calculating confusion matrix
    EvaTest<-predict(model, testing,type="class")
    table(testing$classVar, EvaTest)
    xtab<-table(testing$classVar, EvaTest)
    #Function call
    RPartlist<-ConMat(xtab)
    Results <- rbind(Results, (unlist(RPartlist)))
    #       names(testing)
    ## calculating ROC curves
    #plot(rp)
    q<-predict(model,  newdata=testing,type = "prob")
    pred <- prediction(q[,2], testing$classVar)
    str(pred)
    x[[i]]<-performance(pred, "tpr", "fpr")
    aucc<-performance(pred, "auc")
    aucValues<-rbind(aucValues,aucc@y.values[[1]])
    
    
    #ROC curve
    # png(filename="curve.png")
    
    if (i==1){
      plot(x[[i]], col=i)
    }
    else{
      plot(x[[i]], col=i,add=T)
      abline(0, 1, lty = 2)
      
    }
    #    dev.off()
    
    
  }
  
  
  colnames(aucValues)<-c("AUC")
  Results2<-cbind(Results,aucValues)
  Results3<-round(rbind(Results2,apply(Results2,2,mean)),3)
  colnames(Results3) <- c("Accuracy", "Precision", "Recall", "Specificity", "F-measure","AUC")
  rownames(Results3) <- c(1:N_fold2, "mean")
  #  FileName<-paste(Pathway,"csv", sep=".")
  write.csv(Results3, file="res.CSV")
  return(Results3)
  
}# end of raprt function

############################
#########################

RF_model<-function(N_fold,AGP_data, SedVal){
  # AGP_data<-NewData
  N_fold2<-N_fold
  outcomeName<-names(AGP_data[2])
  
  x = list() # to store ROC curve
  aucValues<-c() # to save auc value
  Results<-data.frame () # to save results from confusion
  set.seed(SedVal)
  folds<-createFolds(AGP_data[,outcomeName], k = N_fold2, list = TRUE, returnTrain = FALSE)
  
  for (i in 1:N_fold2){
    
    # split data into training and testing chunks
    
    testing<- AGP_data[folds[[i]],]
    dim(testing)
    #exclude lowe confidence genes
    testing<-testing[!(testing$ID  %in% Rem),]
    write.csv(testing,"159.CSV")
    
        testing<-testing[,-1]
    ################
    
    training<- AGP_data[-folds[[i]],]
    training<-training[,-1]
   # write.csv(training,"159_train.CSV")    
    
    ###########
    #model desiging
    model<-randomForest(classVar~., data=training)
    #fancyRpartPlot(model)
    #rpart.plot(model)
    
    #Evaluation
    # Calculating confusion matrix
    EvaTest<-predict(model, testing,type="class")
    table(testing$classVar, EvaTest)
    xtab<-table(testing$classVar, EvaTest)
    #Function call
    RPartlist<-ConMat(xtab)
    Results <- rbind(Results, (unlist(RPartlist)))
    #       names(testing)
    ## calculating ROC curves
    #plot(rp)
    q<-predict(model,  newdata=testing,type = "prob")
    pred <- prediction(q[,2], testing$classVar)
    str(pred)
    x[[i]]<-performance(pred, "tpr", "fpr")
    aucc<-performance(pred, "auc")
    aucValues<-rbind(aucValues,aucc@y.values[[1]])
    
    
    
    
  }

  posi<-c(0.35, 0.27,0.22,0.17,0.12, .07)
  
    
for (k in 1:N_fold2){
  #ROC curve
  # png(filename="curve.png")
  
  if (k==1){
    plot(x[[k]], col=k, lwd=2)
    text(1,posi[k],labels=paste("AUC = ",round(aucValues[k],2),sep=""),adj=1,col=k+1)
  }
  else{
    plot(x[[k]], lwd=2,col=k,add=T)
    text(1,posi[k],labels=paste("AUC = ",round(aucValues[k],2),sep=""),adj=1,col=k+1)
    abline(0, 1, lty = 2)
    
  }
  #    dev.off()
  
}
    
  
  colnames(aucValues)<-c("AUC")
  Results2<-cbind(Results,aucValues)
  Results3<-round(rbind(Results2,apply(Results2,2,mean)),2)
  text(1,posi[6],labels=paste("Mean AUC = ",round(Results3[nrow(Results3),ncol(Results3)],2),sep=""),adj=1,col=1)
  
    colnames(Results3) <- c("Accuracy", "Precision", "Recall", "Specificity", "F-measure","AUC")
  rownames(Results3) <- c(1:N_fold2, "mean")
  #  FileName<-paste(Pathway,"csv", sep=".")
  write.csv(Results3, file="res.CSV")
  return(Results3)
  
}# end of raprt function

################
#   Random forest feature selection function
##################################

# Function to calculate confusion matrix
########################################
ConMat<-function(x){
  # Function created to calculate precision, recall and others stuff
  #Its uses confusionMatrix fucntion from CART package
  #take confusion table as input
  #returen a list which should be unlist before writting to data frame
  
  result<-confusionMatrix(x, mode="everything", positive="ASD") #
  Accuarcy<- result$overall['Accuracy']
  sensi<-result$byClass['Sensitivity']
  speci<-result$byClass['Specificity']
  
  preci <- result$byClass['Precision']    
  reca <- result$byClass['Recall']
  f_measure<-result$byClass['F1']
  #  f_measure <- 2 * ((preci * reca) / (preci + reca))
  str(result)
  
  ff<-list ( Accuarcy, preci, reca, speci, f_measure)
  return (ff)
} # end of confusion matrix function
##########################
##################################
Navie.bayes_model<-function(N_fold,AGP_data, seedVal){
  
  
  #AGP_data<-data
  N_fold2<-N_fold
  outcomeName<-names(AGP_data[2])
  
  x = list() # to store ROC curve
  aucValues<-c() # to save auc value
  Results<-data.frame () # to save results from confusion
  set.seed(seedVal)
  folds<-createFolds(AGP_data[,outcomeName], k = N_fold2, list = TRUE, returnTrain = FALSE)
  
  for (i in 1:N_fold2){
    
    # split data into training and testing chunks
    
    testing<- AGP_data[folds[[i]],]
    dim(testing)
    testing<-testing[!(testing$ID  %in% Rem),]
    #write.csv(testing,"159.CSV")
    
    testing<-testing[,-1]
    
    
    ################
    ##########################
    training<- AGP_data[-folds[[i]],]
    dim(training)
    training<-training[,-1]
    
    ###############################
    ##############################################
    #####################################
    # klaR NaiveBayes
    set.seed(3456)
    # ?NaiveBayes
    # nb_model <- NaiveBayes(training[,-1], training$classVar, usekernel=F)
    ####################################
    nb_model<-naiveBayes(training[,-1], training$classVar, usekernel=F)
    ###########
    #?naiveBayes
    #Evaluation
    # Calculating confusion matrix
    EvaTest<-predict(nb_model, testing, type="class")
    # class(EvaTest); str(EvaTest)
    table(testing$classVar,EvaTest)
    xtab<-table(testing$classVar, EvaTest)
    #Function call
    RPartlist<-ConMat(xtab)
    Results <- rbind(Results, (unlist(RPartlist)))
    
    #Auc value
    q<-predict(nb_model,  newdata=testing,type = "raw")
    str(q)
    pred <- prediction(q[,2], testing$classVar)
    str(pred)
    x[[i]]<-performance(pred, "tpr", "fpr")
    aucc<-performance(pred, "auc")
    aucValues<-rbind(aucValues,aucc@y.values[[1]])
    
    
    
    
  }
  
  
  
  Results2<-cbind(Results,aucValues)
  Results3<-round(rbind(Results2,apply(Results2,2,mean)),3)
  colnames(Results3) <- c("Accuracy", "Precision", "Recall", "Specificity", "F-measure","AUC")
  rownames(Results3) <- c(1:N_fold2, "mean")
  #  FileName<-paste(Pathway,"csv", sep=".")
  write.csv(Results3, file="res.CSV")
  return(Results3)
  
}# end of Naive bayes function

################################################



################################################
svm_model2<-function(N_fold,AGP_data,sdVal){
  #AGP_data<-data
  N_fold2<-N_fold
  outcomeName<-names(AGP_data[1])
  
  x = list() # to store ROC curve
  aucValues<-c() # to save auc value
  Results<-data.frame () # to save results from confusion
  set.seed(sdVal)
  folds<-createFolds(AGP_data[,outcomeName], k = N_fold2, list = TRUE, returnTrain = FALSE)
  
  for (i in 1:N_fold2){
    
    # split data into training and testing chunks
    
    testing<- AGP_data[folds[[i]],]
    dim(testing)
    
    #exclude lowe confidence genes
    testing<-testing[!(testing$ID  %in% Rem),]
  #  write.csv(testing,"159.CSV")
    
    testing<-testing[,-1]
    
    ################
    
    training<- AGP_data[-folds[[i]],]
    training<-training[,-1]
    
    ###########?svm
    #model desiging
    model<-svm(classVar~., data=training, method="C-classification",
               kernel="radial", probability = T, gamma=0.02, epsilon = 0.001)
    #Evaluation
    # Calculating confusion matrix
    
    EvaTest<-predict(model,  newdata=testing,type="class")
    table(testing$classVar, EvaTest)
    xtab<-table(testing$classVar, EvaTest)
    #Function call
    RPartlist<-ConMat(xtab)
    Results <- rbind(Results, (unlist(RPartlist)))
    
    #       names(testing)
    ## calculating ROC curves
    #plot(rp)
    q<-predict(model, testing[,-1],probability = T)
    svmmodel.probs<-attr(q,"probabilities")[,2]
    pred <- prediction(svmmodel.probs, testing[,1])
    
    x[[i]]<-performance(pred, "tpr", "fpr")
    aucc<-performance(pred, "auc")
    aucValues<-rbind(aucValues,aucc@y.values[[1]])
    
    
    #ROC curve
    # png(filename="curve.png")
    
    if (i==1){
      plot(x[[i]], col=i)
    }
    else{
      plot(x[[i]], col=i,add=T)
      abline(0, 1, lty = 2)
      
    }
    #    dev.off()
    
    
  }
  
  
  Results2<-cbind(Results,aucValues)
  Results3<-round(rbind(Results2,apply(Results2,2,mean)),3)
  colnames(Results3) <- c("Accuracy", "Precision", "Recall", "Specificity", "F-measure","AUC")
  rownames(Results3) <- c(1:N_fold2, "mean")
  #  FileName<-paste(Pathway,"csv", sep=".")
  write.csv(Results3, file="res.CSV")
  return(Results3)
  
  
}# end of svm function
#################################
##################################################
linear_svm_model<-function(N_fold,AGP_data, Sval){
  #AGP_data<-data
  N_fold2<-N_fold
  outcomeName<-names(AGP_data[1])
  
  x = list() # to store ROC curve
  aucValues<-c() # to save auc value
  Results<-data.frame () # to save results from confusion
  set.seed(Sval)
  folds<-createFolds(AGP_data[,outcomeName], k = N_fold2, list = TRUE, returnTrain = FALSE)
  
  for (i in 1:N_fold2){
    
    # split data into training and testing chunks
    
    testing<- AGP_data[folds[[i]],]
    dim(testing)
    #exclude lowe confidence genes
    testing<-testing[!(testing$ID  %in% Rem),]
    #  write.csv(testing,"159.CSV")
    
    testing<-testing[,-1]
    
    
    ################
    
    training<- AGP_data[-folds[[i]],]
    training<-training[,-1]
    
    
    ###########
    #model desiging
    linear_model<-svm(classVar~., data=training, method="C-classification",
                      kernel="linear", probability = T, epsilon = 0.001)
    #Evaluation
    # Calculating confusion matrix
    
    EvaTest<-predict(linear_model,  newdata=testing,type="class")
    table(testing$classVar, EvaTest)
    xtab<-table(testing$classVar, EvaTest)
    #Function call
    RPartlist<-ConMat(xtab)
    Results <- rbind(Results, (unlist(RPartlist)))
    
    #       names(testing)
    ## calculating ROC curves
    #plot(rp)
    q<-predict(linear_model, testing[,-1],probability = T)
    svmmodel.probs<-attr(q,"probabilities")[,2]
    pred <- prediction(svmmodel.probs, testing$classVar)
    
    x[[i]]<-performance(pred, "tpr", "fpr")
    aucc<-performance(pred, "auc")
    aucValues<-rbind(aucValues,aucc@y.values[[1]])
    
    
    #ROC curve
    # png(filename="curve.png")
    
    if (i==1){
      plot(x[[i]], col=i)
    }
    else{
      plot(x[[i]], col=i,add=T)
      abline(0, 1, lty = 2)
      
    }
    #    dev.off()
    
    
  }
  
  
  Results2<-cbind(Results,aucValues)
  Results3<-round(rbind(Results2,apply(Results2,2,mean)),3)
  colnames(Results3) <- c("Accuracy", "Precision", "Recall", "Specificity", "F-measure","AUC")
  rownames(Results3) <- c(1:N_fold2, "mean")
  #  FileName<-paste(Pathway,"csv", sep=".")
  write.csv(Results3, file="res.CSV")
  return(Results3)
  
  
}# end of svm function

#######################################################
#####################################################
###########################################





















