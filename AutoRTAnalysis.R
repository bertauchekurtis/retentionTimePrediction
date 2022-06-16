# autoRT analysis
# kurtis bertauche

predict <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/AutoRT Results/predict/predict/test_evaluate.csv")
predict2 <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/AutoRT Results/predict2/predict2/test_evaluate.csv")
predict3 <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/AutoRT Results/predict3/predict3/test_evaluate.csv")
predict4 <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/AutoRT Results/predict4/predict4/test_evaluate.csv")

predict <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/AutoRT_data2/run1_results/test_evaluate.csv")
predict2 <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/AutoRT_data2/run2_results/test_evaluate.csv")
predict3 <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/AutoRT_data2/run3_results/test_evaluate.csv")
predict4 <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/RScripts/Updated/AutoRT_data2/run4_results/test_evaluate.csv")

predictResiduals <- predict$y - predict$y_pred
predict2Residuals <- predict2$y - predict2$y_pred
predict3Residuals <- predict3$y - predict3$y_pred
predict4Residuals <- predict4$y - predict4$y_pred

# predict analysis
print("RMSE:")
sqrt(mean(((predictResiduals))^2))

print("MAE:")# MAE calculation
mean(abs(predictResiduals))

# calculate 95% error window size
q <- quantile(predictResiduals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(predict$y, predict$y_pred)


# predict2 analysis
print("RMSE:")
sqrt(mean(((predict2Residuals))^2))

print("MAE:")# MAE calculation
mean(abs(predict2Residuals))

# calculate 95% error window size
q <- quantile(predict2Residuals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(predict2$y, predict2$y_pred)


# predict3 analysis
print("RMSE:")
sqrt(mean(((predict3Residuals))^2))

print("MAE:")# MAE calculation
mean(abs(predict3Residuals))

# calculate 95% error window size
q <- quantile(predict3Residuals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(predict3$y, predict3$y_pred)


# predict4 analysis
print("RMSE:")
sqrt(mean(((predict4Residuals))^2))

print("MAE:")# MAE calculation
mean(abs(predict4Residuals))

# calculate 95% error window size
q <- quantile(predict4Residuals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(predict4$y, predict4$y_pred)