# results for DEEP LC
# kurtis bertauche

#withCalibration <- read.csv(file = "C:/Users/Kurtis/Downloads/PEPREC_DEEPLC_PREDICTIONS_CALIBRATED.csv")
#noCalibration <- read.csv(file = "C:/Users/Kurtis/Downloads/PEPREC_DEEPLC_PREDICTIONS_NO_CALIBRATE.csv")
#useTR <- read.csv(file = "C:/Users/Kurtis/Downloads/PEPREC_DEEPLC_PREDICTIONS_USETR.csv")

withCalibration <- read.csv(file = "C:/Users/Kurtis/Downloads/PEPREC_testing_DATA_TWO_deeplc_predictions_calibration.csv")
useTR <- read.csv(file = "C:/Users/Kurtis/Downloads/PEPREC_testing_DATA_TWO_deeplc_predictions_useTR.csv")

# analysis for with calibration
# RMSE calculation
withCalResiduals <- withCalibration$tr - withCalibration$predicted_tr

print("RMSE:")
sqrt(mean(((withCalResiduals))^2))

print("MAE:")# MAE calculation
mean(abs(withCalResiduals))

# calculate 95% error window size
q <- quantile(withCalResiduals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(withCalibration$tr, withCalibration$predicted_tr)

# analysis for no calibration
# RMSE calculation
noCalResiduals <- noCalibration$tr - noCalibration$predicted_tr

print("RMSE: ")
sqrt(mean(((noCalResiduals)) ^ 2))

print("MAE:")
mean(abs(noCalResiduals))

# calculate 95% error window size
q <- quantile(noCalResiduals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(noCalibration$tr, noCalibration$predicted_tr)

# analysis for use tr
# RMSE calculation
TRResiduals <- useTR$tr - useTR$predicted_tr

print("RMSE: ")
sqrt(mean(((TRResiduals)) ^ 2))

print("MAE:")
mean(abs(TRResiduals))

# calculate 95% error window size
q <- quantile(TRResiduals, probs =c(.025,.975))
abs(q[1]) + abs(q[2]) # total length of window

# correlation
cor(useTR$tr, useTR$predicted_tr)

