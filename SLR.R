# SLR 
# kurtis bertauche
# updated 9 june 2022

dataOne_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_ONE.csv")
dataOne_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_ONE.csv")
dataTwo_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_TWO.csv")
dataTwo_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_TWO.csv")

set.seed(37)

# make models
slr_one <- lm(RetentionTime 
              ~unmodA+unmodC+unmodD+unmodE+unmodF
              +unmodG+unmodH+unmodI+unmodK+unmodL
              +unmodM+unmodN+unmodP+unmodQ+unmodR
              +unmodS+unmodT+unmodV+unmodW+unmodY
              +modS+modT+modY+modM, data = dataOne_train)
slr_two <- lm(RetentionTime 
              ~unmodA+unmodC+unmodD+unmodE+unmodF
              +unmodG+unmodH+unmodI+unmodK+unmodL
              +unmodM+unmodN+unmodP+unmodQ+unmodR
              +unmodS+unmodT+unmodV+unmodW+unmodY
              +modS+modT+modY+modM, data = dataTwo_train)

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

# test model one with data two
calcStats(dataOne_test$RetentionTime, predict(slr_one, dataOne_test))
# test model two with data one
calcStats(dataTwo_test$RetentionTime, predict(slr_two, dataTwo_test))
# get p-values
summary(slr_one)
summary(slr_two)