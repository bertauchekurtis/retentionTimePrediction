# kurtis bertauche
# extreme gradient boosting

data <- read.csv(file = "C:/Users/kbertauche/Downloads/RetentionTime_HCD_Marx2013_SuppT3.csv")
data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")

data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/dataSetTwoFiltered.csv")
data$X <- NULL
colnames(data) <- c("Peptide.Sequence2", "RetentionTime")

set.seed(37) 
library(stringr)
library(xgboost)

# predictor variables
data$peptideLength <- nchar(data$Peptide.Sequence2)

data$unmodA <- str_count(data$Peptide.Sequence2, "A")
data$unmodC <- str_count(data$Peptide.Sequence2, "C")
data$unmodD <- str_count(data$Peptide.Sequence2, "D")
data$unmodE <- str_count(data$Peptide.Sequence2, "E")
data$unmodF <- str_count(data$Peptide.Sequence2, "F")

data$unmodG <- str_count(data$Peptide.Sequence2, "G")
data$unmodH <- str_count(data$Peptide.Sequence2, "H")
data$unmodI <- str_count(data$Peptide.Sequence2, "I")
data$unmodK <- str_count(data$Peptide.Sequence2, "K")
data$unmodL <- str_count(data$Peptide.Sequence2, "L")

data$unmodM <- str_count(data$Peptide.Sequence2, "M")
data$unmodN <- str_count(data$Peptide.Sequence2, "N")
data$unmodP <- str_count(data$Peptide.Sequence2, "P")
data$unmodQ <- str_count(data$Peptide.Sequence2, "Q")
data$unmodR <- str_count(data$Peptide.Sequence2, "R")

data$unmodS <- str_count(data$Peptide.Sequence2, "S")
data$unmodT <- str_count(data$Peptide.Sequence2, "T")
data$unmodV <- str_count(data$Peptide.Sequence2, "V")
data$unmodW <- str_count(data$Peptide.Sequence2, "W")
data$unmodY <- str_count(data$Peptide.Sequence2, "Y")

data$modS <- str_count(data$Peptide.Sequence2, "s")
data$modT <- str_count(data$Peptide.Sequence2, "t")
data$modY <- str_count(data$Peptide.Sequence2, "y")
data$modM <- str_count(data$Peptide.Sequence2, "m")

# remove non-numeric data
data$Peptide.Sequence2 <- NULL

# split data into training/testing sets
setAssignments <- sample(1:2, size = nrow(data), prob = c(0.8, 0.2), replace = TRUE)
trainingData <- data[setAssignments == 1,]
testingData <- data[setAssignments == 2,]

# exract retention times (output) from data
trainingRetentionTimesLabels <- trainingData$RetentionTime
trainingData$RetentionTime <- NULL

# xg boost matrix creation
trainingxgMatrix <- xgb.DMatrix(data.matrix(trainingData), label = trainingRetentionTimesLabels)

# for testing
testingDataLabels <- testingData$RetentionTime
testingData$RetentionTime <- NULL
testingxgMatrix <- xgb.DMatrix(data.matrix(testingData), label = testingDataLabels)

# file to save results
fileLabelsDF <- data.frame(gamma = numeric(),
                           child_weight = numeric(),
                           depth = numeric(),
                           subsample = numeric(),
                           col_subsample = numeric(),
                           eta = numeric(),
                           n_iterations = numeric(),
                           cv_rmse = numeric(),
                           std_deviation = numeric())
#write.csv(fileLabelsDF, 
 #  file = "C:/Users/kbertauche/Downloads/ExtremeGradientBoostingResults.csv",
 #  append = TRUE)

# a matrix to hold hyperparameter combinations
matrixToTry <- matrix(,nrow=0,ncol=6)
for (gamma in c(0, 0.1, 0.2, 0.3, 0.4))
{
  for (child_weight in c(1,2,3,4,5,6))
  {
    for (col_subsample in c(0.8, 0.9, 1))
    {
      for (max_depth in c(9, 10, 11))
      {
        for (subsample in c(0.8, 0.9, 1))
          {
            for (eta in c(0.01, 0.05, 0.08))
            {
      matrixToTry <- rbind(matrixToTry,
                           c(gamma,child_weight,
                             max_depth, subsample, col_subsample,
                             eta))
            }
        }
      }
    }
  }
}

# tune the models
for (row in 1:nrow(matrixToTry))
{
  set.seed(37)
  model <- xgb.cv(booster = "gbtree",
                  objective = "reg:squarederror",
                  gamma = matrixToTry[row, 1],
                  child_weight = matrixToTry[row, 2],
                  max_depth = matrixToTry[row, 3],
                  subsample = matrixToTry[row, 4],
                  col_subsample = matrixToTry[row, 5],
                  eta = matrixToTry[row, 6],
                  nrounds = 10000,
                  nthreads = 28,
                  nfold = 5,
                  print_every_n = 2500,
                  early_stopping_rounds = 2,
                  data = trainingxgMatrix)
  newResult = data.frame(matrixToTry[row, 1],
                         matrixToTry[row, 2],
                         matrixToTry[row, 3],
                         matrixToTry[row, 4],
                         matrixToTry[row, 5],
                         matrixToTry[row, 6],
                         model$best_iteration,
                         min(model$evaluation_log$test_rmse_mean),
                         model$evaluation_log$test_rmse_std[model$best_iteration])
  names(newResult) <- c("gamma",
                        "child_weight",
                        "depth",
                        "subsample",
                        "col_subsample",
                        "eta",
                        "n_iterations",
                        "cv_rmse",
                        "std_deviation")
  write.table(newResult,
              file = "C:/Users/kbertauche/Downloads/ExtremeGradientBoostingResults.csv",
              append = TRUE,
              col.names = FALSE,
              sep = ",")
  rm(model)
}

# analyze things
set.seed(37)
bestModel <-xgboost(booster = "gbtree",
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
                   data = trainingxgMatrix)

predictions <- predict(bestModel, testingxgMatrix)
residual <- testingDataLabels - predictions
print("RMSE:") # RMSE calculation
sqrt(mean(((residual))^2))

print("MAE:")# MAE calculation
mean(abs(residual))

# calculate 95% error window size
q <- quantile(residual, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

cor(testingDataLabels, predictions)
