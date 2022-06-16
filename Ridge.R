# ridge cross
# kurtis bertauche
# updated 9 june 2022

dataOne_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_ONE.csv")
dataOne_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_ONE.csv")
dataTwo_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_TWO.csv")
dataTwo_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_TWO.csv")

set.seed(37)
library(caret)
library(stringr)
library(stats)
library(glmnet)

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

foldid_one <- sample(rep(seq(5), length.out = nrow(dataOne_train)))
foldid_two <- sample(rep(seq(5), length.out = nrow(dataTwo_train)))

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

fit_ridge_cv_mse_one <- cv.glmnet(Xtrain_one, dataOne_train$RetentionTime, alpha = 0, nfolds = 5, foldid = foldid_one)
fit_ridge_cv_mse_two <- cv.glmnet(Xtrain_two, dataTwo_train$RetentionTime, alpha = 0, nfolds = 5, foldid = foldid_two)

# model one
bestModelminLambda_one <- glmnet(x = Xtrain_one,
                                 y = dataOne_train$RetentionTime,
                                 lambda = fit_ridge_cv_mse_one$lambda.min,
                                 alpha = 0)
bestModelOneSELambda_one <- glmnet(x = Xtrain_one,
                                   y = dataOne_train$RetentionTime,
                                   lambda = fit_ridge_cv_mse_one$lambda.1se,
                                   alpha = 0)

# model two
bestModelminLambda_two <- glmnet(x = Xtrain_two,
                                 y = dataTwo_train$RetentionTime,
                                 lambda = fit_ridge_cv_mse_two$lambda.min,
                                 alpha = 0)
bestModelOneSELambda_two <- glmnet(x = Xtrain_two,
                                   y = dataTwo_train$RetentionTime,
                                   lambda = fit_ridge_cv_mse_two$lambda.1se,
                                   alpha = 0)

# testing for min lambda - data one
Xtest_one = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                           unmodG+unmodH+unmodI+unmodK+unmodL+
                           unmodM+unmodN+unmodP+unmodQ+unmodR+
                           unmodS+unmodT+unmodV+unmodW+unmodY+
                           modS+modY+modT+modM+peptideLength, dataOne_test)[, -1]
predictions_one <- predict(bestModelminLambda_one, newx = Xtest_one)
calcStats(dataOne_test$RetentionTime, predictions_one)
# 1se lambda - data one
predictions_one <- predict(bestModelOneSELambda_one, newx = Xtest_one)
calcStats(dataOne_test$RetentionTime, predictions_one)

# testing for min lambda - data two
Xtest_two = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                           unmodG+unmodH+unmodI+unmodK+unmodL+
                           unmodM+unmodN+unmodP+unmodQ+unmodR+
                           unmodS+unmodT+unmodV+unmodW+unmodY+
                           modS+modY+modT+modM+peptideLength, dataTwo_test)[, -1]
predictions_two <- predict(bestModelminLambda_two, newx = Xtest_two)
calcStats(dataTwo_test$RetentionTime, predictions_two)
# 1se lambda - data two
predictions_two <- predict(bestModelOneSELambda_two, newx = Xtest_two)
calcStats(dataTwo_test$RetentionTime, predictions_two)