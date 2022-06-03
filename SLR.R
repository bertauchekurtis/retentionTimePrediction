# kurtis bertauche
# 2 7 2021
# Simple Linear Regression

# data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")
data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/dataSetTwoFiltered.csv")
data$X <- NULL
colnames(data) <- c("Peptide.Sequence2", "RetentionTime")
set.seed(37) 
library(caret)
library(stringr)

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

# slr
lm1 <- lm(RetentionTime 
          ~unmodA+unmodC+unmodD+unmodE+unmodF
          +unmodG+unmodH+unmodI+unmodK+unmodL
          +unmodM+unmodN+unmodP+unmodQ+unmodR
          +unmodS+unmodT+unmodV+unmodW+unmodY
          +modS+modT+modY+modM, data = trainingData)

# analysis
summary(lm1)
residual <- resid(lm1) 
print("RMSE:") # RMSE calculation
sqrt(mean(((testingData$RetentionTime - predict(lm1, newdata = testingData))^2)))

print("MAE:")# MAE calculation
mean(abs(testingData$RetentionTime - predict(lm1, newdata = testingData)))

# calculate 95% error window size
residuals <- testingData$RetentionTime - predict(lm1, newdata = testingData)
q <- quantile(residuals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(predict(lm1, newdata = testingData), testingData$RetentionTime)
