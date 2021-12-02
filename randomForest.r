# random forest
# kurtis bertauche

# load libraries
library(stringr)
library(randomForest)
library(ranger)

# read data set
data <- read.csv(file = "C:/Users/kbertauche/Downloads/RetentionTime_HCD_Marx2013_SuppT3.csv")

# create info about data
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
data$peptideLength <- nchar(data$Peptide.Sequence2)

data$Peptide.Sequence2 <- NULL

set.seed(37)

# split data
setAssignments <- sample(1:2, size = nrow(data), prob = c(0.8, 0.2), replace = TRUE)
trainingData <- data[setAssignments == 1,]
testingData <- data[setAssignments == 2,]

# a matrix to store hyperparameter combos
matrixToTry <- matrix(,nrow=0,ncol=2)

# fil the matrix
for (numTrees in c(500, 750, 1000,2000, 3000, 5000, 10000))
{
  for (mtry in c(6, 8, 10, 12, 14, 16))
  {
    matrixToTry <- rbind(matrixToTry, c(numTrees, mtry))
  }
}

# somewhere to store results
resultdf <- data.frame(numTrees = numeric(),
                       mtry = numeric(),
                       OOB_mse = numeric())

fileLabelsDF <- data.frame(numtrees = numeric(),
                           mtry = numeric(),
                           OOB_mse = numeric())

write.csv(fileLabelsDF, 
          file = "C:/Users/kbertauche/Downloads/resultsRFFixed.csv",
          append = TRUE)

for(row in 1:nrow(matrixToTry))
{
  set.seed(37)
  rfModel <- ranger(
    formula   = RetentionTime ~ ., 
    data      = trainingData, 
    num.trees = matrixToTry[row, 1],
    mtry      = matrixToTry[row, 2],
    min.node.size = 5,
    num.threads = 24
  )
  print(rfModel)
  newResult <- data.frame(
    matrixToTry[row, 1],
    matrixToTry[row, 2],
    rfModel$prediction.error
  )
  rm(rfModel)
  names(newResult) <- c("numTrees",
                        "mtry",
                        "OOB_mse")
  write.table(newResult, 
              file = "C:/Users/kbertauche/Downloads/resultsRFFixed.csv",
              append = TRUE,
              col.names = FALSE,
              sep = ",")
  resultdf <- rbind(resultdf, newResult)
}

# find RMSE of best model
set.seed(37)
bestModel <- ranger(
  formula   = RetentionTime ~ ., 
  data      = trainingData, 
  num.trees = 10000,
  mtry      = 14,
  min.node.size = 5,
  num.threads = 28)
bestModel
predictions <- predict(bestModel, testingData)
differences <- predictions$predictions - testingData$RetentionTime
diffSquare <- (differences * differences)
summ <- sum(diffSquare)
summ <- summ / 18090
sqrt(summ) # final RMSE value
