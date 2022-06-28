# kurtis bertauche
# random forest

dataOne_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_ONE.csv")
dataOne_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_ONE.csv")
dataTwo_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_TWO.csv")
dataTwo_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_TWO.csv")

set.seed(37) 
library(randomForest)
library(stringr)
library(ranger)

dataOne_test$Peptide.Sequence2 <- NULL
dataOne_train$Peptide.Sequence2 <- NULL
dataTwo_test$Peptide.Sequence2 <- NULL
dataTwo_train$Peptide.Sequence2 <- NULL

# create an empty matrix that will hold the parameter combinations for tuning
matrixToTry <- matrix(,nrow=0,ncol=2)

# fill the parameter matrix with the combinations that will be tried
for (numTrees in c(500, 750, 1000,2000, 3000, 5000, 10000))
{
  for (mtry in c(6, 8, 10, 12, 14, 16))
  {
    matrixToTry <- rbind(matrixToTry, c(numTrees, mtry))
  }
}

# creating a data frame to store results of tuning
resultdf_one <- data.frame(numTrees = numeric(),
                       mtry = numeric(),
                       OOB_mse = numeric())

# creating a file to write results to
fileLabelsDF_one <- data.frame(numtrees = numeric(),
                           mtry = numeric(),
                           OOB_mse = numeric())

# creating a data frame to store results of tuning
resultdf_two <- data.frame(numTrees = numeric(),
                           mtry = numeric(),
                           OOB_mse = numeric())

# creating a file to write results to
fileLabelsDF_two <- data.frame(numtrees = numeric(),
                               mtry = numeric(),
                               OOB_mse = numeric())

# write the column names to the file
write.csv(fileLabelsDF_one, 
          file = "C:/Users/kbertauche/Downloads/resultsRFFixed.csv",
          append = TRUE)

# write the column names to the file (data 2)
write.csv(fileLabelsDF_two,
          file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/tuning2/rfResultsData2",
          append = TRUE)

# tune first model
for(row in 1:nrow(matrixToTry))
{
  # for consistency, keep seed
  set.seed(37)
  # rf model call
  rfModel <- ranger(
    formula   = RetentionTime ~ ., 
    data      = dataOne_train, 
    num.trees = matrixToTry[row, 1],
    mtry      = matrixToTry[row, 2],
    min.node.size = 5,
    num.threads = 24
  )
  print(rfModel)
  # store model stats in data frame
  newResult <- data.frame(
    matrixToTry[row, 1],
    matrixToTry[row, 2],
    rfModel$prediction.error
  )
  # clear model from memory
  rm(rfModel)
  # write most recent model's stats to file
  names(newResult) <- c("numTrees",
                        "mtry",
                        "OOB_mse")
  # write to file
  write.table(newResult,
              file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/rfResultsData1.csv",
              append = TRUE,
              col.names = FALSE,
              sep = ",")
  
  # also keep most recent model's stats in enviornment
  resultdf_one <- rbind(resultdf_one, newResult)
}
# model two
for(row in 1:nrow(matrixToTry))
{
  # for consistency, keep seed
  set.seed(37)
  # rf model call
  rfModel <- ranger(
    formula   = RetentionTime ~ ., 
    data      = dataTwo_train, 
    num.trees = matrixToTry[row, 1],
    mtry      = matrixToTry[row, 2],
    min.node.size = 5,
    num.threads = 24
  )
  print(rfModel)
  # store model stats in data frame
  newResult <- data.frame(
    matrixToTry[row, 1],
    matrixToTry[row, 2],
    rfModel$prediction.error
  )
  # clear model from memory
  rm(rfModel)
  # write most recent model's stats to file
  names(newResult) <- c("numTrees",
                        "mtry",
                        "OOB_mse")
  # for data 2
  write.table(newResult,
              file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/rfResultsData2.csv",
              append = TRUE,
              col.names = FALSE,
              sep = ",")
  
  # also keep most recent model's stats in enviornment
  resultdf_two <- rbind(resultdf_two, newResult)
}

set.seed(37)
# rebuild the model
bestModel_one <- ranger(
  formula   = RetentionTime ~ ., 
  data      = dataOne_train, 
  num.trees = 10000,
  mtry      = 14,
  min.node.size = 5,
  num.threads = 28)
bestModel_one

set.seed(37)
# rebuild the model
bestModel_two <- ranger(
  formula   = RetentionTime ~ ., 
  data      = dataTwo_train, 
  num.trees = 5000,
  mtry      = 14,
  min.node.size = 5,
  num.threads = 28)
bestModel_two

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

predictions <- predict(bestModel_one, dataOne_test)
calcStats(dataOne_test$RetentionTime, predictions)

predictions <- predict(bestModel_one, dataTwo_test)
calcStats(dataTwo_test$RetentionTime, predictions)

predictions <- predict(bestModel_two, dataTwo_test)
calcStats(dataTwo_test$RetentionTime, predictions)

predictions <- predict(bestModel_two, dataOne_test)
calcStats(dataOne_test$RetentionTime, predictions)