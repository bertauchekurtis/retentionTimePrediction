# kurtis bertauche
# random forest


# data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")
data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/dataSetTwoFiltered.csv")
data$X <- NULL
colnames(data) <- c("Peptide.Sequence2", "RetentionTime")
set.seed(37) 
library(randomForest)
library(stringr)
library(ranger)

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

# peptide sequence no longer needed
data$Peptide.Sequence2 <- NULL

# split data into trainng/testing sets
setAssignments <- sample(1:2, size = nrow(data), prob = c(0.8, 0.2), replace = TRUE)
trainingData <- data[setAssignments == 1,]
testingData <- data[setAssignments == 2,]

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
resultdf <- data.frame(numTrees = numeric(),
                       mtry = numeric(),
                       OOB_mse = numeric())

# creating a file to write results to
fileLabelsDF <- data.frame(numtrees = numeric(),
                           mtry = numeric(),
                           OOB_mse = numeric())

# write the column names to the file
write.csv(fileLabelsDF, 
          file = "C:/Users/kbertauche/Downloads/resultsRFFixed.csv",
          append = TRUE)

# write the column names to the file (data 2)
write.csv(fileLabelsDF,
          file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/tuning2/rfResultsData2",
          append = TRUE)

# tune the model using OOB MSE on 80% of the dataset - keep 20% behind for testing (same as SLR split)
for(row in 1:nrow(matrixToTry))
{
  # for consistency, keep seed
  set.seed(37)
  
  # rf model call
  rfModel <- ranger(
    formula   = RetentionTime ~ ., 
    data      = trainingData, 
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
  #write.table(newResult, 
  #            file = "C:/Users/kbertauche/Downloads/resultsRFFixed.csv",
  #            append = TRUE,
  #            col.names = FALSE,
  #            sep = ",")
  
  # for data 2
  write.table(newResult,
              file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/rfResultsData2.csv",
              append = TRUE,
              col.names = FALSE,
              sep = ",")
  
  # also keep most recent model's stats in enviornment
  resultdf <- rbind(resultdf, newResult)
}

# find stats of best model

# keep same seed
set.seed(37)
# rebuild the model
bestModel <- ranger(
  formula   = RetentionTime ~ ., 
  data      = trainingData, 
  num.trees = 5000,
  mtry      = 14,
  min.node.size = 5,
  num.threads = 28)
bestModel

# calculate RMSE
predictions <- predict(bestModel, testingData)
residuals <- testingData$RetentionTime - predictions$predictions 


# analysis
print("RMSE:") # RMSE calculation
sqrt(mean(((residuals))^2))

print("MAE:")# MAE calculation
mean(abs(residuals))

# calculate 95% error window size
q <- quantile(residuals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(predictions$predictions, testingData$RetentionTime)

