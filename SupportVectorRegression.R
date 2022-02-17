# kurtis bertauche
# support vector machine

data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")
set.seed(37) 
library(parallel)
library(stringr)
library(e1071)
library(caret)
library(kernlab)

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

# search space for cost parameter
valuesToTune = c(0.1, 0.5, 1, 02, 5, 10, 20)

# function to build model (to be used in parallel)
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
                 data = trainingData, 
                 method = "svmLinear", 
                 trControl = train_control,  
                 tuneGrid=data.frame(C=c(x)))
  
  # collect results
  results <- data.frame("results" = c(x, model$results$RMSE))
  write.csv(results, 
            file = paste("C:/Users/kbertauche/Downloads/SVMrun", x, ".csv",
                         sep="",
                         collapse = NULL))
  model
}

# run model in parallel
cl <- makeCluster(7)
clusterExport(cl, list( 
  "trainControl",
  "trainingData",
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
bestModel <- train(RetentionTime ~
                     unmodA + unmodC + unmodD + unmodE + unmodF +
                     unmodG + unmodH + unmodI + unmodK + unmodL +
                     unmodM + unmodN + unmodP + unmodQ + unmodR +
                     unmodS + unmodT + unmodV + unmodW + unmodY +
                     modS + modT + modY + modM + peptideLength, 
                   data = data, 
                   method = "svmLinear", 
                   C = 0.01)