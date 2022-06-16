# stepwise cross
# kurtis bertauche
# updated 9 june 2022

dataOne_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_ONE.csv")
dataOne_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_ONE.csv")
dataTwo_test <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/testingSet_withVars_DATA_TWO.csv")
dataTwo_train <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/data/trainingSet_withVars_DATA_TWO.csv")

set.seed(37)
library(caret)
library(leaps)

stepwiseModel_one <- regsubsets(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                                  unmodG+unmodH+unmodI+unmodK+unmodL+
                                  unmodM+unmodN+unmodP+unmodQ+unmodR+
                                  unmodS+unmodT+unmodV+unmodW+unmodY+
                                  modS+modY+modT+modM
                                , data = dataOne_train, nvmax = 24, method = "exhaustive")

stepwiseModel_two <- regsubsets(RetentionTime ~unmodA+unmodC+unmodD+unmodE+unmodF+
                                  unmodG+unmodH+unmodI+unmodK+unmodL+
                                  unmodM+unmodN+unmodP+unmodQ+unmodR+
                                  unmodS+unmodT+unmodV+unmodW+unmodY+
                                  modS+modY+modT+modM
                                , data = dataTwo_train, nvmax = 24, method = "exhaustive")


stepwiseModelSum_one <- summary(stepwiseModel_one)
stepwiseModelSum_two <- summary(stepwiseModel_two)

num_vars_one = ncol(dataOne_train) - 3
trn_idx_one = sample(c(TRUE, FALSE), nrow(dataOne_train), replace = TRUE)
tst_idx_one = (!trn_idx_one)

num_vars_two = ncol(dataTwo_train) - 3
trn_idx_two = sample(c(TRUE, FALSE), nrow(dataTwo_train), replace = TRUE)
tst_idx_two = (!trn_idx_two)

test_mat_one = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                              unmodG+unmodH+unmodI+unmodK+unmodL+
                              unmodM+unmodN+unmodP+unmodQ+unmodR+
                              unmodS+unmodT+unmodV+unmodW+unmodY+
                              modS+modY+modT+modM, dataOne_test)

test_mat_two = model.matrix(RetentionTime ~ unmodA+unmodC+unmodD+unmodE+unmodF+
                              unmodG+unmodH+unmodI+unmodK+unmodL+
                              unmodM+unmodN+unmodP+unmodQ+unmodR+
                              unmodS+unmodT+unmodV+unmodW+unmodY+
                              modS+modY+modT+modM, dataTwo_test)

# FOR MODEL ONE
test_err_rmse_one = rep(0, times = num_vars_one)
mae_one = rep(0, times = num_vars_one)
qs_one = rep(0, times = num_vars_one)
cors_one = rep(0, times = num_vars_one)
for (i in seq_along(test_err_rmse_one)) {
  coefs = coef(stepwiseModel_one, id = i)
  pred = test_mat_one[, names(coefs)] %*% coefs
  test_err_rmse_one[i] <- sqrt(mean((dataOne_test$RetentionTime - pred) ^ 2))
  mae_one[i] <- mean(abs(dataOne_test$RetentionTime - pred))
  q <- quantile((dataOne_test$RetentionTime-pred), probs =c(.025,.975))
  qs_one[i] <- abs(q[1]) + abs(q[2]) # total length of window
  cors_one[i] <- cor(pred, dataOne_test$RetentionTime)
}

# FOR MODEL TWO
test_err_rmse_two = rep(0, times = num_vars_two)
mae_two = rep(0, times = num_vars_two)
qs_two = rep(0, times = num_vars_two)
cors_two = rep(0, times = num_vars_two)
for (i in seq_along(test_err_rmse_two)) {
  coefs = coef(stepwiseModel_two, id = i)
  pred = test_mat_two[, names(coefs)] %*% coefs
  test_err_rmse_two[i] <- sqrt(mean((dataTwo_test$RetentionTime - pred) ^ 2))
  mae_two[i] <- mean(abs(dataTwo_test$RetentionTime - pred))
  q <- quantile((dataTwo_test$RetentionTime-pred), probs =c(.025,.975))
  qs_two[i] <- abs(q[1]) + abs(q[2]) # total length of window
  cors_two[i] <- cor(pred, dataTwo_test$RetentionTime)
}

# select best models using RMSE
min_index_one <- which.min(test_err_rmse_one)
min_index_two <- which.min(test_err_rmse_two)

# best model stats (ONE)
test_err_rmse_one[min_index_one]
mae_one[min_index_one]
qs_one[min_index_one]
cors_one[min_index_one]

# best model stats (TWO)
test_err_rmse_two[min_index_two]
mae_two[min_index_two]
qs_two[min_index_two]
cors_two[min_index_two]

# coefs
coef(stepwiseModel_one, min_index_one)
coef(stepwiseModel_two, min_index_two)

# to get p-values
lmodel_one <- lm(dataOne_train$RetentionTime ~ unmodA + unmodC + unmodD + unmodE
                 + unmodF + unmodH + unmodI + unmodK + unmodL +
                   + unmodM +unmodN +  unmodP + unmodQ + unmodR +unmodS + unmodT + unmodV + unmodW + unmodY +
                   + modS + modY + modT + modM, data = dataOne_train)
summary(lmodel_one)

lmodel2 <- lm(dataTwo_train$RetentionTime ~ unmodA + unmodC + unmodD + unmodE
              + unmodF + unmodG + unmodH + unmodI + unmodK + unmodL +
                + unmodM + unmodN + unmodP +unmodQ + unmodR + unmodS + unmodT + unmodV + unmodW + unmodY +
                + modS + modY + modT + modM , data = dataTwo_train)
summary(lmodel2)
