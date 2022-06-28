# kurtis bertauche
# support vector machine
# updated 21 june 2022

dataOne_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_ONE.csv")
dataOne_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_ONE.csv")
dataTwo_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_TWO.csv")
dataTwo_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_TWO.csv")

dataOne_test$Peptide.Sequence2 <- NULL
dataOne_train$Peptide.Sequence2 <- NULL
dataTwo_test$Peptide.Sequence2 <- NULL
dataTwo_train$Peptide.Sequence2 <- NULL


set.seed(37) 
library(parallel)
library(stringr)
library(e1071)
library(caret)
library(kernlab)

# search space for cost parameter
valuesToTune = c(0.1, 0.5, 1, 02, 5, 10, 20)

# FIRST DATA SET
parallelSVM = function(x)
{
  # set seed for consistency
  set.seed(37)
  
  # set cross validation method
  train_control <- trainControl(method="cv", number=5)
  
  # create the model
  model <- train(RetentionTime ~
                   unmodA + unmodC + unmodD + unmodE + unmodF +
                   unmodG + unmodH + unmodI + unmodK + unmodL +
                   unmodM + unmodN + unmodP + unmodQ + unmodR +
                   unmodS + unmodT + unmodV + unmodW + unmodY +
                   modS + modT + modY + modM + peptideLength, 
                 data = dataOne_train, 
                 method = "svmLinear", 
                 trControl = train_control,  
                 tuneGrid=data.frame(C=c(x)))
  
  # collect results
  results <- data.frame("results" = c(x, model$results$RMSE))
  write.csv(results, 
            file = paste("C:/Users/Kurtis/Desktop/retentionTimePrediction/tuning/SVM/SVMrun", x, ".csv",
                         sep="",
                         collapse = NULL))
  model
}

# run model in parallel
cl <- makeCluster(7)
clusterExport(cl, list( 
  "trainControl",
  "dataOne_train",
  "valuesToTune",
  "parallelSVM",
  "train" 
), 
envir=environment())

finalResults <- parSapply(cl = cl, 
                          X = valuesToTune,
                          FUN = parallelSVM)
stopCluster(cl)

# SECOND DATA SET
parallelSVM = function(x)
{
  # set seed for consistency
  set.seed(37)
  
  # set cross validation method
  train_control <- trainControl(method="cv", number=5)
  
  # create the model
  model <- train(RetentionTime ~
                   unmodA + unmodC + unmodD + unmodE + unmodF +
                   unmodG + unmodH + unmodI + unmodK + unmodL +
                   unmodM + unmodN + unmodP + unmodQ + unmodR +
                   unmodS + unmodT + unmodV + unmodW + unmodY +
                   modS + modT + modY + modM + peptideLength, 
                 data = dataTwo_train, 
                 method = "svmLinear", 
                 trControl = train_control,  
                 tuneGrid=data.frame(C=c(x)))
  
  # collect results
  results <- data.frame("results" = c(x, model$results$RMSE))
  write.csv(results, 
            file = paste("C:/Users/Kurtis/Desktop/retentionTimePrediction/tuning2/SVM/SVMrun", x, ".csv",
                         sep="",
                         collapse = NULL))
  model
}

# run model in parallel
cl <- makeCluster(7)
clusterExport(cl, list( 
  "trainControl",
  "dataTwo_train",
  "valuesToTune",
  "parallelSVM",
  "train" 
), 
envir=environment())

finalResults <- parSapply(cl = cl, 
                          X = valuesToTune,
                          FUN = parallelSVM)
stopCluster(cl)


# analysis (to be run after all models have been built and best model is known)
bestModelOne <- train(RetentionTime ~
                     unmodA + unmodC + unmodD + unmodE + unmodF +
                     unmodG + unmodH + unmodI + unmodK + unmodL +
                     unmodM + unmodN + unmodP + unmodQ + unmodR +
                     unmodS + unmodT + unmodV + unmodW + unmodY +
                     modS + modT + modY + modM + peptideLength, 
                   data = dataOne_train, 
                   method = "svmLinear", 
                   cost = 0.1)

bestModelTwo <- train(RetentionTime ~
                        unmodA + unmodC + unmodD + unmodE + unmodF +
                        unmodG + unmodH + unmodI + unmodK + unmodL +
                        unmodM + unmodN + unmodP + unmodQ + unmodR +
                        unmodS + unmodT + unmodV + unmodW + unmodY +
                        modS + modT + modY + modM + peptideLength, 
                      data = dataTwo_train, 
                      method = "svmLinear", 
                      cost = 0.01)

calcStats = function(trueResponse, predictedResponse)
{
  residuals <- trueResponse - predictedResponse
  # RMSE
  rmse <- sqrt(mean(residuals ^ 2))
  # mae
  mae <- mean(abs(residuals))
  # window
  q <- quantile(residuals, probs =c(.025,.975))
  window <- abs(q[1]) + abs(q[2])
  # correlation
  corr <- cor(predictedResponse, trueResponse)
  # return vector
  c(rmse, mae, window, corr)
}

predicts_m1d1 <- predict(bestModelOne, dataOne_test)
predicts_m1d2 <- predict(bestModelOne, dataTwo_test)
predicts_m2d1 <- predict(bestModelTwo, dataOne_test)
predicts_m2d2 <- predict(bestModelTwo, dataTwo_test)

calcStats(dataOne_test$RetentionTime, predicts_m1d1)
calcStats(dataOne_test$RetentionTime, predicts_m2d1)
calcStats(dataTwo_test$RetentionTime, predicts_m1d2)
calcStats(dataTwo_test$RetentionTime, predicts_m2d2)
