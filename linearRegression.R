# linearRegression.r
# kurtis bertauche

##################
# this script contains the code for all of the linear models

##################


data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")
set.seed(37) 
library(caret)
library(stringr)
library(leaps)
library(stats)
library(glmnet)

# split data into trainng/testing sets
setAssignments <- sample(1:2, size = nrow(data), prob = c(0.8, 0.2), replace = TRUE)
trainingData <- data[setAssignments == 1,]
testingData <- data[setAssignments == 2,]

# make a historam to visualize data
trainingRetentionTimes <- trainingData$RetentionTime
hist(trainingRetentionTimes, main = "Training Data", xlab = "Retention Time", col = "#316D9E")

# predictor variables
trainingData$peptideLength <- nchar(trainingData$Peptide.Sequence2)
testingData$peptideLength <- nchar(testingData$Peptide.Sequence2)

trainingData$unmodA <- str_count(trainingData$Peptide.Sequence2, "A")
trainingData$unmodC <- str_count(trainingData$Peptide.Sequence2, "C")
trainingData$unmodD <- str_count(trainingData$Peptide.Sequence2, "D")
trainingData$unmodE <- str_count(trainingData$Peptide.Sequence2, "E")
trainingData$unmodF <- str_count(trainingData$Peptide.Sequence2, "F")

testingData$unmodA <- str_count(testingData$Peptide.Sequence2, "A")
testingData$unmodC <- str_count(testingData$Peptide.Sequence2, "C")
testingData$unmodD <- str_count(testingData$Peptide.Sequence2, "D")
testingData$unmodE <- str_count(testingData$Peptide.Sequence2, "E")
testingData$unmodF <- str_count(testingData$Peptide.Sequence2, "F")

trainingData$unmodG <- str_count(trainingData$Peptide.Sequence2, "G")
trainingData$unmodH <- str_count(trainingData$Peptide.Sequence2, "H")
trainingData$unmodI <- str_count(trainingData$Peptide.Sequence2, "I")
trainingData$unmodK <- str_count(trainingData$Peptide.Sequence2, "K")
trainingData$unmodL <- str_count(trainingData$Peptide.Sequence2, "L")

testingData$unmodG <- str_count(testingData$Peptide.Sequence2, "G")
testingData$unmodH <- str_count(testingData$Peptide.Sequence2, "H")
testingData$unmodI <- str_count(testingData$Peptide.Sequence2, "I")
testingData$unmodK <- str_count(testingData$Peptide.Sequence2, "K")
testingData$unmodL <- str_count(testingData$Peptide.Sequence2, "L")

trainingData$unmodM <- str_count(trainingData$Peptide.Sequence2, "M")
trainingData$unmodN <- str_count(trainingData$Peptide.Sequence2, "N")
trainingData$unmodP <- str_count(trainingData$Peptide.Sequence2, "P")
trainingData$unmodQ <- str_count(trainingData$Peptide.Sequence2, "Q")
trainingData$unmodR <- str_count(trainingData$Peptide.Sequence2, "R")

testingData$unmodM <- str_count(testingData$Peptide.Sequence2, "M")
testingData$unmodN <- str_count(testingData$Peptide.Sequence2, "N")
testingData$unmodP <- str_count(testingData$Peptide.Sequence2, "P")
testingData$unmodQ <- str_count(testingData$Peptide.Sequence2, "Q")
testingData$unmodR <- str_count(testingData$Peptide.Sequence2, "R")

trainingData$unmodS <- str_count(trainingData$Peptide.Sequence2, "S")
trainingData$unmodT <- str_count(trainingData$Peptide.Sequence2, "T")
trainingData$unmodV <- str_count(trainingData$Peptide.Sequence2, "V")
trainingData$unmodW <- str_count(trainingData$Peptide.Sequence2, "W")
trainingData$unmodY <- str_count(trainingData$Peptide.Sequence2, "Y")

testingData$unmodS <- str_count(testingData$Peptide.Sequence2, "S")
testingData$unmodT <- str_count(testingData$Peptide.Sequence2, "T")
testingData$unmodV <- str_count(testingData$Peptide.Sequence2, "V")
testingData$unmodW <- str_count(testingData$Peptide.Sequence2, "W")
testingData$unmodY <- str_count(testingData$Peptide.Sequence2, "Y")

trainingData$modS <- str_count(trainingData$Peptide.Sequence2, "s")
trainingData$modT <- str_count(trainingData$Peptide.Sequence2, "t")
trainingData$modY <- str_count(trainingData$Peptide.Sequence2, "y")
trainingData$modM <- str_count(trainingData$Peptide.Sequence2, "m")

testingData$modS <- str_count(testingData$Peptide.Sequence2, "s")
testingData$modT <- str_count(testingData$Peptide.Sequence2, "t")
testingData$modY <- str_count(testingData$Peptide.Sequence2, "y")
testingData$modM <- str_count(testingData$Peptide.Sequence2, "m")

# scatterplot of data
plot(trainingData$peptideLength, trainingData$RetentionTime, main = "Peptide Length vs. Retention Time",
     xlab = "Peptide Length", ylab = "Retention Time", col = "#316D9E")

###############
#
# LINEAR MODEL
#
###############
lm1 <- lm(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF
          +unmodG+unmodH+unmodI+unmodK+unmodL
          +unmodM+unmodN+unmodP+unmodQ+unmodR
          +unmodS+unmodT+unmodV+unmodW+unmodY
          +modS+modT+modY+peptideLength+modM, data = trainingData)

# evaluate linear regression
predictions <- predict(lm1, testingData)
differences <- testingData$RetentionTime - predictions
differecnes <- differences * differences
me <- mean(differecnes)
sqrt(me)
summary(lm1)
residual <- resid(lm1) 

# plot the residuals
plot(trainingData$RetentionTime, residual, ylab = "Residuals", xlab = "Retention Time",
     main = "Retention Time Residuals (Training)", col = "#316D9E")
abline(0, 0)

# plot the residuals color-coded for modified
plot(trainingData$RetentionTime, residual, ylab = "Residuals", xlab = "Retention Time",
     main = "Retention Time Residuals (Training)", col = ifelse((trainingData$modS + trainingData$modY + trainingData$modT) > 0, "red", "#316D9E"))
abline(0, 0)

# plot predicted vs actual for training
plot(predict(lm1), trainingData$RetentionTime,
     main = "Predicted vs. Actual (Training)",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#316D9E")
abline(a = 0,                                        
       b = 1,
       lwd = 2)

# plot predicted vs actual for testing
plot(predict(lm1, newdata = testingData), testingData$RetentionTime,
     main = "Predicting vs. Actual (Testing)",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#088572")
abline(a = 0,
       b = 1,
       lwd = 2)

# plot predicted vs actual for training with modified highlighted
plot(predict(lm1), trainingData$RetentionTime,
     main = "Predicted vs. Actual (Training)",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = ifelse((trainingData$modS + trainingData$modY + trainingData$modT) > 0, "red", "#316D9E"))
abline(a = 0,                                        
       b = 1,
       lwd = 2)

# plot predicted vs actual for testing with modified highlighted
plot(predict(lm1, newdata = testingData), testingData$RetentionTime,
     main = "Predicted Values vs. Actual Values (Testing Set)",
     xlab = "Predicted Values (min)",
     ylab = "Observed Values (min)",
     col = ifelse((testingData$modS + testingData$modY + testingData$modT) > 0, "red", "#088572"))
abline(a = 0,
       b = 1,
       lwd = 2)

# manual calculation of RMSE
sqrt(mean(((testingData$RetentionTime - predict(lm1, newdata = testingData))^2)))

###############
#
# STEPWISE REGRESSION
#
###############

# RMSE Function
# inputs: actual values and predicted values
# output: RMSE value for given inputs
rmse = function(actual, predicted) {
        sqrt(mean((actual - predicted) ^ 2))
}

# RMSE Function (more versatile)
# inputs: a model, input data, response variable
# output: RMSE value for a model on the input dataset for the response variable
get_rmse = function(model, data, response) {
        rmse(actual = data[, response], 
             predicted = predict(model, data))
}

# Complexity Function (# of predictor variables)
# inputs: a model
# output: compelxity of a model
get_complexity = function(model) {
        length(coef(model)) - 1
}

# naive stepwise regression
# aribtrarily choosen predictor variables added one by one
lm2 <- lm(RetentionTime ~peptideLength, data = trainingData)
lm3 <- lm(RetentionTime ~unmodA, data = trainingData)
lm4 <- lm(RetentionTime ~unmodA+unmodC, data = trainingData)
lm5 <- lm(RetentionTime ~unmodA+unmodC+unmodD, data = trainingData)
lm6 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE, data = trainingData)
lm7 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF, data = trainingData)
lm8 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                        unmodG, data = trainingData)
lm9 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                  unmodG+unmodH, data = trainingData)
lm10 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                    unmodG+unmodH+unmodI, data = trainingData)
lm11 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK, data = trainingData)
lm12 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL, data = trainingData)
lm13 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                unmodG+unmodH+unmodI+unmodK+unmodL+
                unmodM, data = trainingData)
lm14 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN, data = trainingData)
lm15 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP, data = trainingData)
lm16 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ, data = trainingData)
lm17 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR, data = trainingData)
lm18 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS, data = trainingData)
lm19 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT, data = trainingData)
lm20 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT+unmodV, data = trainingData)
lm21 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT+unmodV+unmodW, data = trainingData)
lm22 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT+unmodV+unmodW+unmodY, data = trainingData)
lm23 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT+unmodV+unmodW+unmodY+
                   modS, data = trainingData)
lm24 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT+unmodV+unmodW+unmodY+
                   modS+modY, data = trainingData)
lm25 <- lm(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                   unmodG+unmodH+unmodI+unmodK+unmodL+
                   unmodM+unmodN+unmodP+unmodQ+unmodR+
                   unmodS+unmodT+unmodV+unmodW+unmodY+
                   modS+modY+modT, data = trainingData)

# analysis of naive stepwise regression
modelList <- list(lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11, lm12,
                  lm13, lm14, lm15, lm16, lm17, lm18, lm19, lm20, lm21, lm22,
                  lm23, lm24, lm25, lm1)
trainRmse = sapply(modelList, get_rmse, data = trainingData, response = "RetentionTime")
testRmse = sapply(modelList, get_rmse, data = testingData, response = "RetentionTime")
model_complexity = sapply(modelList, get_complexity)

# plots model complexity (size) vs. RMSE for training and test sets
plot(model_complexity, trainRmse, type = "b", 
     ylim = c(min(c(trainRmse, testRmse)) - 0.02, 
              max(c(trainRmse, testRmse)) + 0.02), 
     col = "#316D9E", 
     xlab = "Model Size",
     ylab = "RMSE")
lines(model_complexity, testRmse, type = "b", col = "#088572")

# output the RMSE for testing and training sets
testRmse
trainRmse

# formal stepwise regression
fit_all2 = regsubsets(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                             unmodG+unmodH+unmodI+unmodK+unmodL+
                             unmodM+unmodN+unmodP+unmodQ+unmodR+
                             unmodS+unmodT+unmodV+unmodW+unmodY+
                             modS+modY+modT+modM+
                             peptideLength , data = trainingData, nvmax = 24)

# analysis
fitAllSum = summary(fit_all2)
names(fitAllSum)
fitAllSum$bic

# plot the rss and highlight best (lowest)
plot(fitAllSum$rss, xlab = "Number of Variables", ylab = "RSS")
bestRSS = which.min(fitAllSum$rss)
points(bestRSS, fitAllSum$rss[bestRSS], col = "red", cex = 2, pch = 20)

# plot the adjr and highlight the best (highest)
plot(fitAllSum$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq")
bestAdjR2 = which.max(fitAllSum$adjr2)
points(bestAdjR2, fitAllSum$adjr2[bestAdjR2], col = "red", cex = 2, pch = 20)

# plot mallows cp and highlight the best (lowest)
plot(fitAllSum$cp, xlab = "Number of Variables", ylab = "Cp")
bestCP = which.min(fitAllSum$cp)
points(bestCP, fitAllSum$cp[bestCP], col = "red", cex = 2, pch = 20)

# plot bic and highlight the best (lowest)
plot(fitAllSum$bic, xlab = "Number of Variables", ylab = "BIC")
bestBIC = which.min(fitAllSum$bic)
points(bestBIC, fitAllSum$bic[bestBIC], col = "red", cex = 2, pch = 20)

# forward stepwise regression
fitFwd = regsubsets(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                            unmodG+unmodH+unmodI+unmodK+unmodL+
                            unmodM+unmodN+unmodP+unmodQ+unmodR+
                            unmodS+unmodT+unmodV+unmodW+unmodY+modM+
                            modS+modY+modT+peptideLength, 
                    data = trainingData, method = "forward")
fitFwdSum = summary(fitFwd)

# backwards method stepwise regression
fitBwd = regsubsets(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                            unmodG+unmodH+unmodI+unmodK+unmodL+
                            unmodM+unmodN+unmodP+unmodQ+unmodR+
                            unmodS+unmodT+unmodV+unmodW+unmodY+modM+
                            modS+modY+modT+peptideLength, 
                    data = trainingData, method = "backward")

# analysis
fitBwdSum = summary(fitBwd)
which.min(fitBwdSum$cp)
coef(fitBwd, which.min(fitBwdSum$cp))
fitAicBack = step(lm, trace = FALSE)
coef(fitAicBack)

# validated (test set) RMSE
testDataMatrix <- model.matrix(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                                       unmodG+unmodH+unmodI+unmodK+unmodL+
                                       unmodM+unmodN+unmodP+unmodQ+unmodR+
                                       unmodS+unmodT+unmodV+unmodW+unmodY+
                                       modS+modY+modT+modM+peptideLength, data = testingData)

testErr <- rep(0, times = (ncol(testingData)-4))
for(i in seq_along(testErr))
{
        coef(fit_all2, 1)
        i
        coefs <- coef(fit_all2, i)
        pred = testDataMatrix[, names(coefs)] %*% coefs
        testErr[i] <- sqrt(mean((testingData$RetentionTime - pred) ^ 2))
}
testErr
plot(testErr, type='b', ylab = "Test Set RMSE", xlab = "Number of Predictors")
bestTestRMSE = which.min(testErr)
points(bestTestRMSE, testErr[bestTestRMSE], col = "red", cex = 2, pch = 20)
coef(fit_all2, 22)
which.min(testErr)

# cross validated RMSE (adapted from R for Statistical Learning) ####

numFolds <- 5
data2 <- data

data2$unmodA <- str_count(data2$Peptide.Sequence2, "A")
data2$unmodC <- str_count(data2$Peptide.Sequence2, "C")
data2$unmodD <- str_count(data2$Peptide.Sequence2, "D")
data2$unmodE <- str_count(data$Peptide.Sequence2, "E")
data2$unmodF <- str_count(data$Peptide.Sequence2, "F")

data2$unmodG <- str_count(data$Peptide.Sequence2, "G")
data2$unmodH <- str_count(data$Peptide.Sequence2, "H")
data2$unmodI <- str_count(data$Peptide.Sequence2, "I")
data2$unmodK <- str_count(data$Peptide.Sequence2, "K")
data2$unmodL <- str_count(data$Peptide.Sequence2, "L")

data2$unmodM <- str_count(data$Peptide.Sequence2, "M")
data2$unmodN <- str_count(data$Peptide.Sequence2, "N")
data2$unmodP <- str_count(data$Peptide.Sequence2, "P")
data2$unmodQ <- str_count(data$Peptide.Sequence2, "Q")
data2$unmodR <- str_count(data$Peptide.Sequence2, "R")

data2$unmodS <- str_count(data$Peptide.Sequence2, "S")
data2$unmodT <- str_count(data$Peptide.Sequence2, "T")
data2$unmodV <- str_count(data$Peptide.Sequence2, "V")
data2$unmodW <- str_count(data$Peptide.Sequence2, "W")
data2$unmodY <- str_count(data$Peptide.Sequence2, "Y")

data2$modS <- str_count(data$Peptide.Sequence2, "s")
data2$modY <- str_count(data$Peptide.Sequence2, "y")
data2$modT <- str_count(data$Peptide.Sequence2, "t")
data2$modM <- str_count(data$Peptide.Sequence2, "m")
data2$peptideLength <- nchar(data$Peptide.Sequence2)

predict.regsubsets = function(object, newdata, id, ...) {
        
        form  = as.formula(object$call[[2]])
        mat   = model.matrix(form, newdata)
        coefs = coef(object, id = id)
        xvars = names(coefs)
        
        mat[, xvars] %*% coefs
}

folds = caret::createFolds(data2$RetentionTime, k = numFolds)
fold_error = matrix(0, nrow = numFolds, ncol = 23, 
                    dimnames = list(paste(1:5), paste(1:23)))

for(j in 1:numFolds) {
        
        train_fold    = data2[-folds[[j]], ]
        validate_fold = data2[ folds[[j]], ]
        
        best_fit = regsubsets(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                                      unmodG+unmodH+unmodI+unmodK+unmodL+
                                      unmodM+unmodN+unmodP+unmodQ+unmodR+
                                      unmodS+unmodT+unmodV+unmodW+unmodY+
                                      modS+modY+modT+modM+peptideLength, data = train_fold, nvmax = 23)
        for (i in 1:23) {
                pred = predict(best_fit, validate_fold, id = i)
                fold_error[j, i] = rmse(actual = validate_fold$RetentionTime,                                        predicted = pred)
        }
}

cv_error = apply(fold_error, 2, mean)
cv_error

# plots cross-validated RMSE vs. number of predictor variables
plot(cv_error, type='b', ylab = "Corss-Validated RMSE", xlab = "Number of Predictors")
bestCVRMSE = which.min(cv_error)
points(bestCVRMSE, cv_error[bestCVRMSE], col = "red", cex = 2, pch = 20)

# plots predicted vs. actual values for testing set
plot(predict(fit_all2, testingData, id = 23), testingData$RetentionTime,
     main = "Predicted vs. Actual (Testing)",
     sub = "Stepwise Regression",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#088572")
abline(a = 0,
       b = 1,
       lwd = 2)

###############
# ridge regression
# adapted from R for statistcal learning
#
###############

# creating input values for glmnet function
X = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                         unmodG+unmodH+unmodI+unmodK+unmodL+
                         unmodM+unmodN+unmodP+unmodQ+unmodR+
                         unmodS+unmodT+unmodV+unmodW+unmodY+
                         modS+modY+modT+modM+peptideLength, trainingData)[, -1]
Xtest = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                         unmodG+unmodH+unmodI+unmodK+unmodL+
                         unmodM+unmodN+unmodP+unmodQ+unmodR+
                         unmodS+unmodT+unmodV+unmodW+unmodY+
                         modS+modY+modT+modM+peptideLength, testingData)[, -1]
y = trainingData$RetentionTime
lmm <- lm(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                  unmodG+unmodH+unmodI+unmodK+unmodL+
                  unmodM+unmodN+unmodP+unmodQ+unmodR+
                  unmodS+unmodT+unmodV+unmodW+unmodY+
                  modS+modY+modT+modM+peptideLength, data = trainingData)
coef(lmm)
sum(abs(coef(lmm)[-1]))
sum(coef(lmm)[-1] ^ 2)

# fits the ridge regression
fit_ridge = glmnet(X, y, alpha = 0)
# plots L1 Norm
plot(fit_ridge)
# plots Log Lambda
plot(fit_ridge, xvar = "lambda", label = TRUE)

# cross validated ridge regression
fit_ridge_cv = cv.glmnet(X, y, alpha = 0)

# plots log of lambda vs. MSE
plot(fit_ridge_cv)
# extract results
coef(fit_ridge_cv)
coef(fit_ridge_cv, s = "lambda.min")
sum(coef(fit_ridge_cv, s = "lambda.1se")[-1] ^ 2) # penalty term using 1-SE rule lambda

# predict using minimum lambda
predict(fit_ridge_cv, X, s = "lambda.min")

# plots predicted vs actual (training) using lambda min
plot(predict(fit_ridge_cv, X, s = "lambda.min"), trainingData$RetentionTime,
     main = "Predicted vs. Actual (Training)",
     sub = "Ridge Regression - Lambda Min",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#316D9E")
abline(a = 0,
       b = 1,
       lwd = 2)
# plots predicted vs actual (testing) using lambda min
plot(predict(fit_ridge_cv, Xtest, s = "lambda.min"), testingData$RetentionTime,
     main = "Predicted vs. Actual (Testing)",
     sub = "Ridge Regression - Lambda Min",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#088572")
abline(a = 0,
       b = 1,
       lwd = 2)

# predict using lambda se
predict(fit_ridge_cv, X)

# plots predicted vs actual (training) using lambda se
plot(predict(fit_ridge_cv, X), trainingData$RetentionTime,
     main = "Predicted vs. Actual (Training)",
     sub = "Ridge Regression - Lambda SE",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#316D9E")
abline(a = 0,
       b = 1,
       lwd = 2)
# plots predicted vs actual (teting) using lambda se
plot(predict(fit_ridge_cv, Xtest), testingData$RetentionTime,
     main = "Predicted vs. Actual (Testing)",
     sub = "Ridge Regression - Lambda SE",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#088572")
abline(a = 0,
       b = 1,
       lwd = 2)

# calculate MSE
mean((y - predict(fit_ridge_cv, X)) ^ 2)
# calculate RMSE
sqrt(mean((y - predict(fit_ridge_cv, X)) ^ 2)) 

# CV-RMSEs
sqrt(fit_ridge_cv$cvm)

# CV-RMSE using minimum lambda
sqrt(fit_ridge_cv$cvm[fit_ridge_cv$lambda == fit_ridge_cv$lambda.min])

# CV-RMSE using 1-SE rule lambda
sqrt(fit_ridge_cv$cvm[fit_ridge_cv$lambda == fit_ridge_cv$lambda.1se]) 

###############
# lasso regression
# adapted from R for Statistical Learning
#
###############

fit_lasso = glmnet(X, y, alpha = 1)
plot(fit_lasso)
plot(fit_lasso, xvar = "lambda", label = TRUE)

fit_lasso_cv = cv.glmnet(X, y, alpha = 1)
plot(fit_lasso_cv)
coef(fit_lasso_cv)
coef(fit_lasso_cv, s = "lambda.min")

#cross validated rmse with min lambda:
sqrt(fit_lasso_cv$cvm[fit_lasso_cv$lambda == fit_lasso_cv$lambda.min])
sqrt(fit_lasso_cv$cvm[fit_lasso_cv$lambda == fit_lasso_cv$lambda.1se]) 

# plots predicted vs atucal (training) for lamda min
plot(predict(fit_lasso_cv, X, s = "lambda.min"), trainingData$RetentionTime,
     main = "Predicted vs. Actual (Training)",
     sub = "Lasso Regression - Lambda Min",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#316D9E")
abline(a = 0,
       b = 1,
       lwd = 2)
# plots predicted vs actual (testing) for lambda min
plot(predict(fit_lasso_cv, Xtest, s = "lambda.min"), testingData$RetentionTime,
     main = "Predicted vs. Actual (Testing)",
     sub = "Lasso Regression - Lambda Min",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#088572")
abline(a = 0,
       b = 1,
       lwd = 2)

# plos predicted vs actual (training) for lambda se
plot(predict(fit_lasso_cv, X), trainingData$RetentionTime,
     main = "Predicted vs. Actual (Training)",
     sub = "Lasso Regression - Lambda SE",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#316D9E")
abline(a = 0,
       b = 1,
       lwd = 2)
# plots predicted vs actual (testing) for lambda se
plot(predict(fit_lasso_cv, Xtest), testingData$RetentionTime,
     main = "Predicted vs. Actual (Testing)",
     sub = "Lasso Regression - Lambda SE",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#088572")
abline(a = 0,
       b = 1,
       lwd = 2)

###############
# elastic net
# adapted from R for statistical learning
#
###############

cv_5 = trainControl(method = "cv", number = 5)
elasticNet = train(
        RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                unmodG+unmodH+unmodI+unmodK+unmodL+
                unmodM+unmodN+unmodP+unmodQ+unmodR+
                unmodS+unmodT+unmodV+unmodW+unmodY+
                modS+modY+modT+modM+peptideLength, data = data2,
        method = "glmnet",
        trControl = cv_5
)
elasticNet

elasticNetBig = train(
        RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                unmodG+unmodH+unmodI+unmodK+unmodL+
                unmodM+unmodN+unmodP+unmodQ+unmodR+
                unmodS+unmodT+unmodV+unmodW+unmodY+
                modS+modY+modT+modM+peptideLength ^ 2, data = data2,
        method = "glmnet",
        trControl = cv_5,
        tuneLength = 10
)


get_best_result = function(caret_fit) {
        best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
        best_result = caret_fit$results[best, ]
        rownames(best_result) = NULL
        best_result
}
get_best_result(elasticNetBig)

# plots predicted vs actual (training) for lambda se
plot(predict(elasticNet, X), trainingData$RetentionTime,
     main = "Predicted vs. Actual (Training)",
     sub = "Elastic Net",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#316D9E")
abline(a = 0,
       b = 1,
       lwd = 2)
# plots predicted vs actual (testing) for lambda min
plot(predict(elasticNet, Xtest, s = "lambda.min"), testingData$RetentionTime,
     main = "Predicted vs. Actual (Testing)",
     sub = "Elastic Net",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#088572")
abline(a = 0,
       b = 1,
       lwd = 2)
# plots predicted vs actual (training) for lambda se expanded tuning
plot(predict(elasticNetBig, X), trainingData$RetentionTime,
     main = "Predicted vs. Actual (Training)",
     sub = "Elastic Net - Expanded Tuning",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#316D9E")
abline(a = 0,
       b = 1,
       lwd = 2)

# plots predicted vs actual (testing) for lambda se expanded tuning
plot(predict(elasticNetBig, Xtest), testingData$RetentionTime,
     main = "Predicted vs. Actual (Testing)",
     sub = "Elastic Net - Expanded Tuning",
     xlab = "Predicted Values",
     ylab = "Observed Values",
     col = "#088572")
abline(a = 0,
       b = 1,
       lwd = 2)