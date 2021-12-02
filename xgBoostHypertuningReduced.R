# xg boost hyperparameter tuning
# kurtis bertauche
# 11 october 2021

# load libraries
library(stringr)
library(xgboost)

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

# remove non-numeric data 
data$Peptide.Sequence2 <- NULL

# exract retention times (output) from data
retentionTimesLabels <- data$RetentionTime
data$RetentionTime <- NULL

# put data in an xgboost matrix
xgMatrix <- xgb.DMatrix(data.matrix(data), label = retentionTimesLabels)

# make file
fileLabelsDF <- data.frame(n_estimators = numeric(),
                           gamma = numeric(),
                           child_weight = numeric(),
                           depth = numeric(),
                           subsample = numeric(),
                           col_subsample = numeric(),
                           eta = numeric(),
                           n_iterations = numeric(),
                           cv_rmse = numeric(),
                           std_deviation = numeric())

write.csv(fileLabelsDF, 
          file = "C:/Users/kbertauche/Downloads/results4.csv",
          append = TRUE)

ReducedmatrixToTry <- matrix(,nrow=0,ncol=7)
for (n_estimators in c(300, 330, 350, 370, 400))
{
  for (gamma in c(0, 0.1, 0.2, 0.3, 0.4))
  {
    for (child_weight in c(1, 2, 3, 4, 5, 6))
    {
      for (col_subsample in c(0.6, 0.7, 0.8, 0.9))
      {
          ReducedmatrixToTry <- rbind(ReducedmatrixToTry,
                                c(n_estimators,gamma,child_weight,
                                  11, 0.9, col_subsample,
                                  0.05))
      }
    }
  }
}

for (row in 1:nrow(ReducedmatrixToTry))
{
  set.seed(37)
  model <- xgb.cv(booster = "gbtree",
                  objective = "reg:squarederror",
                  n_estimators = ReducedmatrixToTry[row, 1],
                  gamma = ReducedmatrixToTry[row, 2],
                  child_weight = ReducedmatrixToTry[row, 3],
                  max_depth = ReducedmatrixToTry[row, 4],
                  subsample = ReducedmatrixToTry[row, 5],
                  col_subsample = ReducedmatrixToTry[row, 6],
                  eta = ReducedmatrixToTry[row, 7],
                  nrounds = 10000,
                  nthreads = 28,
                  nfold = 5,
                  print_every_n = 2500,
                  early_stopping_rounds = 2,
                  data = xgMatrix)
  newResult <- data.frame(ReducedmatrixToTry[row ,1],
                          ReducedmatrixToTry[row, 2],
                          ReducedmatrixToTry[row, 3],
                          ReducedmatrixToTry[row, 4],
                          ReducedmatrixToTry[row, 5],
                          ReducedmatrixToTry[row, 6],
                          ReducedmatrixToTry[row, 7],
                          model$best_iteration,
                          min(model$evaluation_log$test_rmse_mean),
                          model$evaluation_log$test_rmse_std[model$best_iteration])
  names(newResult) <- c("n_estimators",
                        "gamma",
                        "child_weight",
                        "depth",
                        "subsample",
                        "col_subsample",
                        "eta",
                        "n_iterations",
                        "cv_rmse",
                        "std_deviation")
  write.table(newResult, 
              file = "C:/Users/kbertauche/Downloads/results4.csv",
              append = TRUE,
              col.names = FALSE,
              sep = ",")
  rm(model)
}
                                          