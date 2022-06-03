# kurtis bertauche
# 2 9 2021
# Stepwise Linear Regression

 data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")
data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/dataSetTwoFiltered.csv")
data$X <- NULL
colnames(data) <- c("Peptide.Sequence2", "RetentionTime")
set.seed(37) 
library(caret)
library(stringr)
library(leaps)

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

# split data into trainng/testing sets
setAssignments <- sample(1:2, size = nrow(data), prob = c(0.8, 0.2), replace = TRUE)
trainingData <- data[setAssignments == 1,]
testingData <- data[setAssignments == 2,]

stepwiseModel <- regsubsets(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                        unmodG+unmodH+unmodI+unmodK+unmodL+
                        unmodM+unmodN+unmodP+unmodQ+unmodR+
                        unmodS+unmodT+unmodV+unmodW+unmodY+
                        modS+modY+modT+modM
                         , data = trainingData, nvmax = 24, method = "exhaustive")
stepwiseModelSum <- summary(stepwiseModel)

coef(stepwiseModel, id = 1
     )

num_vars = ncol(trainingData) - 3
trn_idx = sample(c(TRUE, FALSE), nrow(trainingData), replace = TRUE)
tst_idx = (!trn_idx)

test_mat = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                       unmodG+unmodH+unmodI+unmodK+unmodL+
                       unmodM+unmodN+unmodP+unmodQ+unmodR+
                       unmodS+unmodT+unmodV+unmodW+unmodY+
                       modS+modY+modT+modM+peptideLength, testingData)


coe <- coef(stepwiseModel, id = 1)

test_err_rmse = rep(0, times = num_vars)
mae = rep(0, times = num_vars)
qs = rep(0, times = num_vars)
cors = rep(0, times = num_vars)
for (i in seq_along(test_err_rmse)) {
  coefs = coef(stepwiseModel, id = i)
  pred = test_mat[, names(coefs)] %*% coefs
  test_err_rmse[i] <- sqrt(mean((testingData$RetentionTime - pred) ^ 2))
  mae[i] <- mean(abs(testingData$RetentionTime - pred))
  q <- quantile((testingData$RetentionTime-pred), probs =c(.025,.975))
  qs[i] <- abs(q[1]) + abs(q[2]) # total length of window
  cors[i] <- cor(pred, testingData$RetentionTime)
}


# all of the test RMSEs
test_err_rmse
# all of the test MAEs
mae
# all of the windows
qs
# all of the correlations
cors

min_index <- which.min(test_err_rmse) # select for model based on rmse

# best model stats
test_err_rmse[min_index]
mae[min_index]
qs[min_index]
cors[min_index]

coef(stepwiseModel, min_index)

# create the best model to get p-values of coefficeints (data 1)
lmodel <- lm(trainingData$RetentionTime ~ unmodA + unmodC + unmodD + unmodE
            + unmodF + unmodH + unmodI + unmodK + unmodL +
              + unmodM +unmodN +  unmodP + unmodQ + unmodR +unmodS + unmodT + unmodV + unmodW + unmodY +
              + modS + modY + modT + modM, data = trainingData)
summary(lmodel)

# creat the best model to get p-values of coefficeients (data 2)
lmodel2 <- lm(trainingData$RetentionTime ~ unmodA + unmodC + unmodD + unmodE
             + unmodF + unmodG + unmodG + unmodH + unmodI + unmodK + unmodL +
               + unmodM + unmodN + unmodP +unmodQ + unmodR + unmodS + unmodT + unmodV + unmodW + unmodY +
               + modS + modY + modT + modM , data = trainingData)
summary(lmodel2)
