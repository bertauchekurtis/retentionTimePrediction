# align
# 22 june 2022
# kurtis bertauche

set.seed(37)
library(glmnet)
library(ranger)
library(randomForest)
library(xgboost)
library(e1071)
library(kernlab)
library(caret)
dataOne_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_ONE.csv")
dataOne_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_ONE.csv")
dataTwo_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_TWO.csv")
dataTwo_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_TWO.csv")

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

# parameter model A is any model created
# dataB is the data that needs to be aligned to it
alignStats = function(modelA, dataB_test, dataB_train)
{
  # get predictor var x for lm - predicted rts of B's train set using model A
  predicts_of_train_B_from_model_A <- predict(modelA, dataB_train)
  # create df with predictor var x and response var y (actual rts of B's train set)
  df <- data.frame(dataB_train$RetentionTime, predicts_of_train_B_from_model_A)
  # create the LM with previously made df
  newLM <- lm(dataB_train.RetentionTime ~ predicts_of_train_B_from_model_A, data = df)
  
  # now, predict the rts of the B's test set with Model A
  predicts_of_test_B_from_model_A <- predict(modelA, dataB_test)
  # previous line will be the input, so makde into df
  df2 <- data.frame(predicts_of_test_B_from_model_A)
  # chnage in name to avoid conflict with predict function
  colnames(df2) <- c("predicts_of_train_B_from_model_A")
  
  # predict the actual rts of B's test using the predicts of B from model A
  predictions_z <- predict(newLM, df2)
  # return analysis
  c(newLM$coefficients, calcStats(dataB_test$RetentionTime, predictions_z))
}

# Simple Linear Regression
slr_one <- lm(RetentionTime 
              ~unmodA+unmodC+unmodD+unmodE+unmodF
              +unmodG+unmodH+unmodI+unmodK+unmodL
              +unmodM+unmodN+unmodP+unmodQ+unmodR
              +unmodS+unmodT+unmodV+unmodW+unmodY
              +modS+modT+modY+modM, data = dataOne_train)

alignStats(slr_one, dataTwo_test, dataTwo_train)

slr_two <- lm(RetentionTime 
              ~unmodA+unmodC+unmodD+unmodE+unmodF
              +unmodG+unmodH+unmodI+unmodK+unmodL
              +unmodM+unmodN+unmodP+unmodQ+unmodR
              +unmodS+unmodT+unmodV+unmodW+unmodY
              +modS+modT+modY+modM, data = dataTwo_train)
alignStats(slr_two, dataOne_test, dataOne_train)

# Stepwise
bestStepwise_one <- lm(dataOne_train$RetentionTime ~ unmodA + unmodC + unmodD + unmodE
                 + unmodF + unmodH + unmodI + unmodK + unmodL +
                   + unmodM +unmodN +  unmodP + unmodQ + unmodR +unmodS + unmodT + unmodV + unmodW + unmodY +
                   + modS + modY + modT + modM, data = dataOne_train)
alignStats(bestStepwise_one, dataTwo_test, dataTwo_train)

bestStepwise_two <- lm(dataTwo_train$RetentionTime ~ unmodA + unmodC + unmodD + unmodE
              + unmodF + unmodG + unmodH + unmodI + unmodK + unmodL +
                + unmodM + unmodN + unmodP +unmodQ + unmodR + unmodS + unmodT + unmodV + unmodW + unmodY +
                + modS + modY + modT + modM , data = dataTwo_train)
alignStats(bestStepwise_two, dataOne_test, dataOne_train)

# put Data into format for glmnet functions
Xtrain_one <- model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                             unmodG+unmodH+unmodI+unmodK+unmodL+
                             unmodM+unmodN+unmodP+unmodQ+unmodR+
                             unmodS+unmodT+unmodV+unmodW+unmodY+
                             modS+modY+modT+modM+peptideLength, dataOne_train)[, -1]
Xtrain_two <- model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                             unmodG+unmodH+unmodI+unmodK+unmodL+
                             unmodM+unmodN+unmodP+unmodQ+unmodR+
                             unmodS+unmodT+unmodV+unmodW+unmodY+
                             modS+modY+modT+modM+peptideLength, dataTwo_train)[, -1]
Xtest_one = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                           unmodG+unmodH+unmodI+unmodK+unmodL+
                           unmodM+unmodN+unmodP+unmodQ+unmodR+
                           unmodS+unmodT+unmodV+unmodW+unmodY+
                           modS+modY+modT+modM+peptideLength, dataOne_test)[, -1]
Xtest_two = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                           unmodG+unmodH+unmodI+unmodK+unmodL+
                           unmodM+unmodN+unmodP+unmodQ+unmodR+
                           unmodS+unmodT+unmodV+unmodW+unmodY+
                           modS+modY+modT+modM+peptideLength, dataTwo_test)[, -1]
# another function for the data formats that GLMNET expects
alignStatsGLM = function(modelA, dataB_test_matrix, dataB_train_matrix, testRTs, trainRTs)
{
  # get predictor var x for lm - predicted rts of B's train set using model A
  predicts_of_train_B_from_model_A <- predict(modelA, dataB_train_matrix)
  # create df with predictor var x and response var y (actual rts of B's train set)
  df <- data.frame(trainRTs, predicts_of_train_B_from_model_A)
  # create the LM with previously made df
  newLM <- lm(trainRTs ~ s0, data = df)
  
  # now, predict the rts of the B's test set with Model A
  predicts_of_test_B_from_model_A <- predict(modelA, dataB_test_matrix)
  # previous line will be the input, so makde into df
  df2 <- data.frame(predicts_of_test_B_from_model_A)
  # chnage in name to avoid conflict with predict function
  colnames(df2) <- c("s0")
  
  # predict the actual rts of B's test using the predicts of B from model A
  predictions_z <- predict(newLM, df2)
  # return analysis
  c(newLM$coefficients, calcStats(testRTs, predictions_z))
}

# Ridge
# Lambda and Alpha values were retreived from bestHyperparameters.tsv
ridge_minLambda_one <- glmnet(x = Xtrain_one,
                                 y = dataOne_train$RetentionTime,
                                 lambda = 0.5111337,
                                 alpha = 0)
alignStatsGLM(ridge_minLambda_one, Xtest_two, Xtrain_two, dataTwo_test$RetentionTime, dataTwo_train$RetentionTime)
ridge_minLambda_two <- glmnet(x = Xtrain_two,
                                   y = dataTwo_train$RetentionTime,
                                   lambda = 2.25766,
                                   alpha = 0)
alignStatsGLM(ridge_minLambda_two, Xtest_one, Xtrain_one, dataOne_test$RetentionTime, dataOne_train$RetentionTime)
ridge_1se_one <- glmnet(x = Xtrain_one,
                              y = dataOne_train$RetentionTime,
                              lambda = 1.075887,
                              alpha = 0)
alignStatsGLM(ridge_1se_one, Xtest_two, Xtrain_two, dataTwo_test$RetentionTime, dataTwo_train$RetentionTime)
ridge_1se_two <- glmnet(x = Xtrain_two,
                              y = dataTwo_train$RetentionTime,
                              lambda = 3.945325,
                              alpha = 0)
alignStatsGLM(ridge_1se_two, Xtest_one, Xtrain_one, dataOne_test$RetentionTime, dataOne_train$RetentionTime)

# LASSO
lasso_minLambda_one <- glmnet(x = Xtrain_one,
                                 y = dataOne_train$RetentionTime,
                                 lambda = 0.0120857,
                                 alpha = 1)
alignStatsGLM(lasso_minLambda_one, Xtest_two, Xtrain_two, dataTwo_test$RetentionTime, dataTwo_train$RetentionTime)
lasso_minLambda_two <- glmnet(x = Xtrain_two,
                              y = dataTwo_train$RetentionTime,
                              lambda = 0.04038162,
                              alpha = 1)
alignStatsGLM(lasso_minLambda_two, Xtest_one, Xtrain_one, dataOne_test$RetentionTime, dataOne_train$RetentionTime)



lasso_OneSELambda_one <- glmnet(x = Xtrain_one,
                                   y = dataOne_train$RetentionTime,
                                   lambda = 0.1635253,
                                   alpha = 1)
alignStatsGLM(lasso_OneSELambda_one, Xtest_two, Xtrain_two, dataTwo_test$RetentionTime, dataTwo_train$RetentionTime)
lasso_OneSELambda_two <- glmnet(x = Xtrain_two,
                                y = dataTwo_train$RetentionTime,
                                lambda = 0.4133184,
                                alpha = 1)
alignStatsGLM(lasso_OneSELambda_two, Xtest_one, Xtrain_one, dataOne_test$RetentionTime, dataOne_train$RetentionTime)

# ELASTIC NET
elastic_one <- glmnet(x = Xtrain_one,
                        y = dataOne_train$RetentionTime,
                        lambda = 0.01739136,
                        alpha = 0.475)
alignStatsGLM(elastic_one, Xtest_two, Xtrain_two, dataTwo_test$RetentionTime, dataTwo_train$RetentionTime)
elastic_two <- glmnet(x = Xtrain_two,
                        y = dataTwo_train$RetentionTime,
                        lambda = 0.03782381,
                        alpha = 0.8875)
alignStatsGLM(elastic_two, Xtest_one, Xtrain_one, dataOne_test$RetentionTime, dataOne_train$RetentionTime)

# RF
dataOne_test$Peptide.Sequence2 <- NULL
dataOne_train$Peptide.Sequence2 <- NULL
dataTwo_test$Peptide.Sequence2 <- NULL
dataTwo_train$Peptide.Sequence2 <- NULL
set.seed(37)
# rebuild the model
rf_one <- ranger(
  formula   = RetentionTime ~ ., 
  data      = dataOne_train, 
  num.trees = 10000,
  mtry      = 14,
  min.node.size = 5,
  num.threads = 28)


set.seed(37)
# rebuild the model
rf_two <- ranger(
  formula   = RetentionTime ~ ., 
  data      = dataTwo_train, 
  num.trees = 5000,
  mtry      = 14,
  min.node.size = 5,
  num.threads = 28)

alignStatsRF = function(modelA, dataB_test, dataB_train)
{
  # get predictor var x for lm - predicted rts of B's train set using model A
  predicts_of_train_B_from_model_A <- predict(modelA, dataB_train)
  # create df with predictor var x and response var y (actual rts of B's train set)
  df <- data.frame(dataB_train$RetentionTime, predicts_of_train_B_from_model_A$predictions)
  # create the LM with previously made df
  newLM <- lm(dataB_train.RetentionTime ~ predicts_of_train_B_from_model_A.predictions, data = df)
  
  # now, predict the rts of the B's test set with Model A
  predicts_of_test_B_from_model_A <- predict(modelA, dataB_test)
  # previous line will be the input, so makde into df
  df2 <- data.frame(predicts_of_test_B_from_model_A$predictions)
  # chnage in name to avoid conflict with predict function
  colnames(df2) <- c("predicts_of_train_B_from_model_A.predictions")
  
  # predict the actual rts of B's test using the predicts of B from model A
  predictions_z <- predict(newLM, df2)
  # return analysis
  c(newLM$coefficients, calcStats(dataB_test$RetentionTime, predictions_z))
}

alignStats(rf_one, dataTwo_test, dataTwo_train)
alignStats(rf_two, dataOne_test, dataOne_train)

# XGB
# put data into xgbs format
trainingRetentionTimesLabelsOne <- dataOne_train$RetentionTime
trainingRetentionTimesLabelsTwo <- dataTwo_train$RetentionTime

trainingxgMatrixOne <- xgb.DMatrix(data.matrix(dataOne_train[,-1]), label = trainingRetentionTimesLabelsOne)
trainingxgMatrixTwo <- xgb.DMatrix(data.matrix(dataTwo_train[,-1]), label = trainingRetentionTimesLabelsTwo)

testingDataLabelsOne <- dataOne_test$RetentionTime
testingDataLabelsTwo <- dataTwo_test$RetentionTime

testingxgMatrixOne <- xgb.DMatrix(data.matrix(dataOne_test[-1]), label = testingDataLabelsOne)
testingxgMatrixTwo <- xgb.DMatrix(data.matrix(dataTwo_test[,-1]), label = testingDataLabelsTwo)

alignStatsXGB = function(modelA, dataB_test, dataB_train, dataB_test_rts, dataB_train_rts)
{
  # get predictor var x for lm - predicted rts of B's train set using model A
  predicts_of_train_B_from_model_A <- predict(modelA, dataB_train)
  # create df with predictor var x and response var y (actual rts of B's train set)
  df <- data.frame(dataB_train_rts, predicts_of_train_B_from_model_A)
  # create the LM with previously made df
  newLM <- lm(dataB_train_rts ~ predicts_of_train_B_from_model_A, data = df)
  
  # now, predict the rts of the B's test set with Model A
  predicts_of_test_B_from_model_A <- predict(modelA, dataB_test)
  # previous line will be the input, so makde into df
  df2 <- data.frame(predicts_of_test_B_from_model_A)
  # chnage in name to avoid conflict with predict function
  colnames(df2) <- c("predicts_of_train_B_from_model_A")
  
  # predict the actual rts of B's test using the predicts of B from model A
  predictions_z <- predict(newLM, df2)
  # return analysis
  c(newLM$coefficients, calcStats(dataB_test_rts, predictions_z))
}

set.seed(37)
xgb_one <-xgboost(booster = "gbtree",
                    objective = "reg:squarederror",
                    gamma = 0.2,
                    child_weight = 1,
                    max_depth = 11,
                    subsample = 0.9,
                    col_subsample = 1,
                    eta = 0.01,
                    nrounds = 10000,
                    nthreads = 28,
                    print_every_n = 2500,
                    early_stopping_rounds = 2,
                    data = trainingxgMatrixOne)
###############33TEST THIS WITH DATA TWO
predictions <- predict(xgb_one, testingxgMatrixTwo)
calcStats(dataTwo_test$RetentionTime, predictions)

alignStatsXGB(xgb_one, testingxgMatrixTwo, trainingxgMatrixTwo, dataTwo_test$RetentionTime, dataTwo_train$RetentionTime)

# SVM
svmOne <- train(RetentionTime ~
                        unmodA + unmodC + unmodD + unmodE + unmodF +
                        unmodG + unmodH + unmodI + unmodK + unmodL +
                        unmodM + unmodN + unmodP + unmodQ + unmodR +
                        unmodS + unmodT + unmodV + unmodW + unmodY +
                        modS + modT + modY + modM + peptideLength, 
                      data = dataOne_train, 
                      method = "svmLinear", 
                      cost = 0.1)
svmTwo <- train(RetentionTime ~
                        unmodA + unmodC + unmodD + unmodE + unmodF +
                        unmodG + unmodH + unmodI + unmodK + unmodL +
                        unmodM + unmodN + unmodP + unmodQ + unmodR +
                        unmodS + unmodT + unmodV + unmodW + unmodY +
                        modS + modT + modY + modM + peptideLength, 
                      data = dataTwo_train, 
                      method = "svmLinear", 
                      cost = 0.1)

# AUTO RT
Testpredictions_m1d2 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_firstDataModel_testedWithSecondData/test_evaluate.csv")
Testpredictions_m1td2 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_firstDataTransfer_testedWithSecondData/test_evaluate.csv")
Testpredictions_m2d1 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_secondDataModel_testedWithFirstData/test_evaluate.csv")
Testpredictions_m2td2 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_secondDataTransfer_testedwithFirstData/test_evaluate.csv")

Trainpredictions_m1d2 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_TRAINING_firstDataModel_testedWithSecondData/test_evaluate.csv")
Trainpredictions_m1td2 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_TRAINING_firstDataTransfer_testedWithSecondData/test_evaluate.csv")
Trainpredictions_m2d1 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_TRAINING_secondDataModel_testedWithFirstData/test_evaluate.csv")
Trainpredictions_m2td2 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_TRAINING_secondDataTransfer_testedwithFirstData/test_evaluate.csv")

alignAutoRT = function(predicts_of_train_B_from_A, predicts_of_test_B_from_A, actual_test_rts_B, actual_train_rts_B)
{
  # get predictor var x from lm - predicted rt's of B's train set using model A
  predicts_of_train_B_from_model_A <- predicts_of_train_B_from_A
  # create df with predictor var x and response var y (actual rts of B's train set)
  df <- data.frame(actual_train_rts_B, predicts_of_train_B_from_model_A)
  # creat the LM with previously made df
  newLM <- lm(actual_train_rts_B ~ predicts_of_train_B_from_model_A, data = df)
  
  # now, predicts the rts of the B's test set with model A (already done here)
  predicts_of_test_B_from_model_A <- predicts_of_test_B_from_A
  # make into df
  df2 <- data.frame(predicts_of_test_B_from_model_A)
  # change name
  colnames(df2) <- c("predicts_of_train_B_from_model_A")
  
  # predict the acutal RTs of B's test using the predicts of B from model A
  predictions_z <- predict(newLM, df2)
  # return analysis
  c(newLM$coefficients, calcStats(actual_test_rts_B, predictions_z))
}

alignAutoRT(Trainpredictions_m1d2$y_pred, Testpredictions_m1d2$y_pred, Testpredictions_m1d2$y, Trainpredictions_m1d2$y)
alignAutoRT(Trainpredictions_m1td2$y_pred, Testpredictions_m1td2$y_pred, Testpredictions_m1td2$y, Trainpredictions_m1td2$y)
alignAutoRT(Trainpredictions_m2d1$y_pred, Testpredictions_m2d1$y_pred, Testpredictions_m2d1$y, Trainpredictions_m2d1$y)
alignAutoRT(Trainpredictions_m2td2$y_pred, Testpredictions_m2td2$y_pred, Testpredictions_m2td2$y, Trainpredictions_m2td2$y)

# ARD
ardTestPredictions_m2d1 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/ARD/ARD_train2_testPredict1.csv", sep = ",")
ardTrainPredictions_m2d1 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/ARD/ARD_train2_trainPredict1.csv")

ardTestPredictions_m1d2 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/ARD/ARD_train1_testPredict2.csv")
ardTrainPredictions_m1d2 <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/ARD/ARD_train1_trainPredict2.csv")

alignAutoRT(ardTrainPredictions_m2d1$X..y_pred, ardTestPredictions_m2d1$X..y_pred, dataOne_test$RetentionTime, dataOne_train$RetentionTime)
alignAutoRT(ardTrainPredictions_m1d2$X..y_pred, ardTestPredictions_m1d2$X..y_pred, dataTwo_test$RetentionTime, dataTwo_train$RetentionTime)
