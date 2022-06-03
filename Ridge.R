# kurtis bertauche
# 2 9 2022
# ridge regression

# data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")
data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/dataSetTwoFiltered.csv")
data$X <- NULL
colnames(data) <- c("Peptide.Sequence2", "RetentionTime")
set.seed(37) 
library(caret)
library(stringr)
library(stats)
library(glmnet)

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
foldid <- sample(rep(seq(5), length.out = nrow(trainingData)))

# create training data without retention time
Xtrain <- model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT+unmodV+unmodW+unmodY+
                   modS+modY+modT+modM+peptideLength, trainingData)[, -1]

fit_ridge_cv_mse <- cv.glmnet(Xtrain, trainingData$RetentionTime, alpha = 0, nfolds = 5, foldid = foldid)

bestModelminLambda <- glmnet(x = Xtrain,
                             y = trainingData$RetentionTime,
                             lambda = fit_ridge_cv_mse$lambda.min,
                             alpha = 0)

bestModelOneSELambda <- glmnet(x = Xtrain,
                               y = trainingData$RetentionTime,
                               lambda = fit_ridge_cv_mse$lambda.1se,
                               alpha = 0)


# testing for min lambda
Xtest = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT+unmodV+unmodW+unmodY+
                   modS+modY+modT+modM+peptideLength, testingData)[, -1]
predictions <- predict(bestModelminLambda, newx = Xtest)
residuals <- testingData$RetentionTime - predictions

# analysis for min lambda
print("RMSE:") # RMSE calculation
sqrt(mean(((residuals))^2))

print("MAE:")# MAE calculation
mean(abs(residuals))

# calculate 95% error window size
q <- quantile(residuals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(predictions, testingData$RetentionTime)

# testing for lambda 1 SE
predictions <- predict(bestModelOneSELambda, newx = Xtest)
residuals <- testingData$RetentionTime - predictions

# analysis for lambda 1 SE
print("RMSE:") # RMSE calculation
sqrt(mean(((residuals))^2))

print("MAE:")# MAE calculation
mean(abs(residuals))

# calculate 95% error window size
q <- quantile(residuals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(predictions, testingData$RetentionTime)
