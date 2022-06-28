# elastic net
# 16 june 2022
# kurtis bertauche

dataOne_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_ONE.csv")
dataOne_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_ONE.csv")
dataTwo_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_TWO.csv")
dataTwo_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_TWO.csv")

set.seed(37) 
library(caret)
library(stringr)
library(stats)
library(glmnet)

foldid_one <- sample(rep(seq(5), length.out = nrow(dataOne_train)))
foldid_two <- sample(rep(seq(5), length.out = nrow(dataTwo_train)))

cv_5 <- trainControl(method = "cv", number = 5)

# create models
elasticNet_one = train(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                         unmodG+unmodH+unmodI+unmodK+unmodL+
                         unmodM+unmodN+unmodP+unmodQ+unmodR+
                         unmodS+unmodT+unmodV+unmodW+unmodY+
                         modS+modY+modT+modM+peptideLength ^ 2, 
                       data = dataOne_train,
                       method = "glmnet",
                       trControl = cv_5,
                       foldid = foldid_one,
                       tuneLength = 25
)

elasticNet_two = train(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                         unmodG+unmodH+unmodI+unmodK+unmodL+
                         unmodM+unmodN+unmodP+unmodQ+unmodR+
                         unmodS+unmodT+unmodV+unmodW+unmodY+
                         modS+modY+modT+modM+peptideLength ^ 2, 
                       data = dataTwo_train,
                       method = "glmnet",
                       trControl = cv_5,
                       foldid = foldid_two,
                       tuneLength = 25
)
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

bestModel_one <- glmnet(x = Xtrain_one,
                        y = dataOne_train$RetentionTime,
                        lambda = elasticNet_one$bestTune$lambda,
                        alpha = elasticNet_one$bestTune$alpha)

bestModel_two <- glmnet(x = Xtrain_two,
                        y = dataTwo_train$RetentionTime,
                        lambda = elasticNet_two$bestTune$lambda,
                        alpha = elasticNet_two$bestTune$alpha)

Xtest_one <- model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                            unmodG+unmodH+unmodI+unmodK+unmodL+
                            unmodM+unmodN+unmodP+unmodQ+unmodR+
                            unmodS+unmodT+unmodV+unmodW+unmodY+
                            modS+modY+modT+modM+peptideLength, dataOne_test)[, -1]
Xtest_two <- model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                            unmodG+unmodH+unmodI+unmodK+unmodL+
                            unmodM+unmodN+unmodP+unmodQ+unmodR+
                            unmodS+unmodT+unmodV+unmodW+unmodY+
                            modS+modY+modT+modM+peptideLength, dataTwo_test)[, -1]

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

predictions <- predict(bestModel_one, newx = Xtest_one)
calcStats(dataOne_test$RetentionTime, predictions)

predictions <- predict(bestModel_two, newx = Xtest_two)
calcStats(dataTwo_test$RetentionTime, predictions)


# predict data two using model one
predictions <- predict(bestModel_one, newx = Xtest_two)
calcStats(dataTwo_test$RetentionTime, predictions)

# predict data one using model two
predictions <- predict(bestModel_two, newx = Xtest_one)
calcStats(dataOne_test$RetentionTime, predictions)