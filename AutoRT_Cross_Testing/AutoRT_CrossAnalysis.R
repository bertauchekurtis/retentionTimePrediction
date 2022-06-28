# autoRT cross analysis
# 21 june 2022
# kurtis bertauche

modelOne_testedwTwo <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_firstDataModel_testedWithSecondData/test_evaluate.csv")
modelTwo_testedwOne <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_secondDataModel_testedWithFirstData/test_evaluate.csv")
modelOneTransfer_testedwTwo <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_firstDataTransfer_testedWithSecondData/test_evaluate.csv")
modelTwoTransfer_testedwOne <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/AutoRT_Cross_Testing/result_secondDataTransfer_testedwithFirstData/test_evaluate.csv")

calcStats = function(trueResponse, predictedResponse)
{
  residuals <- trueResponse - predictedResponse
  # RMSE
  rmse <- sqrt(mean(residuals ^ 2))
  # mae
  mae <- mean(abs(residuals))
  # window
  q <- quantile(residuals, probs =c(.025,.975))
  window <- abs(q[1]) + abs(q[2])
  # correlation
  corr <- cor(predictedResponse, trueResponse)
  # return vector
  c(rmse, mae, window, corr)
}

# model One tested with Two
calcStats(modelOne_testedwTwo$y, modelOne_testedwTwo$y_pred)
# model two tested with One
calcStats(modelTwo_testedwOne$y, modelTwo_testedwOne$y_pred)
# model one transfer tested with two
calcStats(modelOneTransfer_testedwTwo$y, modelOneTransfer_testedwTwo$y_pred)
# model two transfer tested with one
calcStats(modelTwoTransfer_testedwOne$y, modelTwoTransfer_testedwOne$y_pred)
