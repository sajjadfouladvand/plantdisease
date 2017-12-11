normalization<-function(data_temp, n_method="minmax"){
  if(n_method=="minmax")
    for(i in 1:ncol(data_temp)){
      data_temp[,i]<-(data_temp[,i]-min(data_temp[,i]))/(max(data_temp[,i]-min(data_temp[,i])))
    }
  #else if(method=="zscore")
  return( data_temp) 
}

plot_roc<-function(train_set_all, test_set, test_set_label, lower_bound_cost,change_cost_step , change_cost_ieteration, features, label){
  require(ggplot2)
  require(reshape2)
  cost_val<-lower_bound_cost
  performance_roc<-matrix(0,change_cost_ieteration,4)
  for(l in 1:change_cost_ieteration){
    predictData <- svm(Inoculated ~ ., data = train_set_all, cost=cost_val , type="C", kernel="sigmoid" )
    par.pred <- predict(predictData, test_set)
    true_positive<-0
    true_negative<-0
    false_positive<-0
    false_negative<-0
    accuracy<-0
    sensitivity <-0
    false_p_rate<-0
    for(j in 1:length(test_set_label)){
      if(test_set_label[j]=="No" && par.pred[j]== "No")
        true_negative<-true_negative+1
      else if(test_set_label[j]=="No" && par.pred[j]== "Yes")
        false_positive<-false_positive+1
      else if(test_set_label[j]=="Yes" && par.pred[j]== "Yes")
        true_positive<-true_positive+1
      else if(test_set_label[j]=="Yes" && par.pred[j]== "No")
        false_negative<-false_negative+1
    }
    accuracy<- (true_negative+true_positive)/(true_negative+true_positive+false_negative+false_positive)
    sensitivity<- (true_positive)/(true_positive+false_negative)
    false_p_rate<- (false_positive)/(false_positive+true_negative)
    performance_roc[l,1]<-sensitivity
    performance_roc[l,2]<-false_p_rate
    performance_roc[l,3]<-accuracy
    performance_roc[l,4]<-cost_val
    cost_val<-cost_val+change_cost_step
  }
  performance_roc<-as.data.frame(performance_roc)
  colnames(performance_roc)<-c("sensitivity", "false alarm rate", "accuracy","cost")
  write.csv(performance_roc, "res.csv")
  performance_roc <- melt(performance_roc ,  id.vars = "cost", variable.name = 'Measures')
  ggplot(performance_roc, aes(cost,value)) + geom_line(aes(colour = Measures))
}
main<-function(filename="RN2016.txt",cross_validation_folds=10, lower_bound_cost=0.0001,change_cost_step=0.0001, change_cost_ieteration=500,x = "PAR",y = c("Phi2", "qL", "PhiNPQ", "PhiNO", "LEF", "FmPrime", "Fs", "FoPrime", "SPAD"),label="Inoculated"){
  data <- read.table("C:\\Courses\\Bioinformatics\\Final Project\\RN2016.txt", header = T, sep="")
  data<- na.omit(data)
  #SIMPLIFY IT TO A TWO CLASS CLASSIFICATIONNNNNNNNNNN
  #SIMPLIFY IT TO A TWO CLASS CLASSIFICATIONNNNNNNNNNN
  #SIMPLIFY IT TO A TWO CLASS CLASSIFICATIONNNNNNNNNNN
  
  #cross_validation_folds<-10
  #lower_bound_cost <- 0.0001
  #change_cost_step <- 0.0001
  #change_cost_ieteration <- 500
  
  #x = "PAR"
  #y = c("Phi2", "qL", "PhiNPQ", "PhiNO", "LEF", "FmPrime", "Fs", "FoPrime", "SPAD")
  features <- c(x, y)
  label<-"Inoculated"
  
  #data_label<-data[[label]]
  data_temp<-subset(data, select = c(features))
  data_temp<-normalization(data_temp, n_method="minmax")
  sds<-cbind(data_temp, data[[label]])
  colnames(sds)[ncol(sds)]<- label
  #sds <- subset(data, select = c(features, label))
  
  
  
  library(e1071)
  
  splited_sds<-split(sds, sample(rep(1:3)))
  test_set<-splited_sds$`1`
  train_set_all<-rbind(splited_sds$`2`, splited_sds$`3`)
  #set.seed(10)
  #cv_sds<-split(train_set_all, sample(rep(1:cross_validation_folds)))
  cv_index<-sample(1:cross_validation_folds,nrow(train_set_all),replace = T)
  train_set_all<-cbind(train_set_all, cv_index)
  best_costs<-matrix(0,cross_validation_folds,1)
  for(i in 1:cross_validation_folds){
    cat("Training and evaluating using Fold ",i," .....")
    train_set<-train_set_all[which(train_set_all[,ncol(train_set_all)]!=i),]
    validation_set<-train_set_all[which(train_set_all[,ncol(train_set_all)]==i),]
    cost_val<- lower_bound_cost
    performance<-matrix(0,change_cost_ieteration,4) # first column is sensitivity, second is false positive rate, third is accuracy and the last is cost 
    train_set<- subset(train_set, select = c(x, y, label))
    validation_label<-validation_set[[label]]
    validation_set<-subset(validation_set, select = c(x, y))
    for(l in 1:change_cost_ieteration){
      #print(l)
      predictData <- svm(Inoculated ~ ., data = train_set, cost=cost_val , type="C", kernel="sigmoid" )
      par.pred <- predict(predictData, validation_set)
      true_positive<-0
      true_negative<-0
      false_positive<-0
      false_negative<-0
      accuracy<-0
      sensitivity <-0
      false_p_rate<-0
      for(j in 1:length(validation_label)){
        if(validation_label[j]=="No" && par.pred[j]== "No")
          true_negative<-true_negative+1
        else if(validation_label[j]=="No" && par.pred[j]== "Yes")
          false_positive<-false_positive+1
        else if(validation_label[j]=="Yes" && par.pred[j]== "Yes")
          true_positive<-true_positive+1
        else if(validation_label[j]=="Yes" && par.pred[j]== "No")
          false_negative<-false_negative+1
      }
      accuracy<- (true_negative+true_positive)/(true_negative+true_positive+false_negative+false_positive)
      sensitivity<- (true_positive)/(true_positive+false_negative)
      false_p_rate<- (false_positive)/(false_positive+true_negative)
      performance[l,1]<-sensitivity
      performance[l,2]<-false_p_rate
      performance[l,3]<-accuracy
      performance[l,4]<-cost_val
      cost_val<-cost_val+change_cost_step
    }
    best_costs[i,1]<-mean(performance[which(performance[,3]==max(performance[,3])),4])
  }
  optimized_cost<-mean(best_costs)
  print("Optimization finished! The cost values is:\n")
  print(optimized_cost)
  
  #***********************************
  #***********************************
  #*******TESTING PHASE***************
  #***********************************
  #***********************************
  train_set_all<- subset(train_set_all, select = c(x, y, label))
  test_set_label<-test_set[[label]]
  test_set<-subset(test_set, select = c(x, y))
  
  
  predictData <- svm(Inoculated ~ ., data = train_set_all, cost=optimized_cost , type="C", kernel="sigmoid" )
  par.pred <- predict(predictData, test_set)
  true_positive<-0
  true_negative<-0
  false_positive<-0
  false_negative<-0
  accuracy<-0
  sensitivity <-0
  false_p_rate<-0
  for(j in 1:length(test_set_label)){
    if(test_set_label[j]=="No" && par.pred[j]== "No")
      true_negative<-true_negative+1
    else if(test_set_label[j]=="No" && par.pred[j]== "Yes")
      false_positive<-false_positive+1
    else if(test_set_label[j]=="Yes" && par.pred[j]== "Yes")
      true_positive<-true_positive+1
    else if(test_set_label[j]=="Yes" && par.pred[j]== "No")
      false_negative<-false_negative+1
  }
  accuracy<- (true_negative+true_positive)/(true_negative+true_positive+false_negative+false_positive)
  sensitivity<- (true_positive)/(true_positive+false_negative)
  false_p_rate<- (false_positive)/(false_positive+true_negative)
  print("*******************")
  print("*******************")
  print("*******************")
  print("*******************")
  print("Testing is finished. Testing performance is:")
  cat("Accuracy is:", accuracy)
  print("\n")
  cat("sensitivity is: ", sensitivity)
  print("\n")
  cat("False positive rate is: ", false_p_rate)
  print("\n")
  plot_roc(train_set_all, test_set, test_set_label, lower_bound_cost,change_cost_step , change_cost_ieteration, features, label)
}
