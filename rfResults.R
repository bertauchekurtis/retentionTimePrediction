# kurtis bertauche
# rf visualization

# load data
data <- read.csv(file = "C:/Users/Kurtis/Desktop/retentionTimePrediction/rfResultsData2.csv", header =FALSE)
data$V1 <- NULL
colnames(data) <- c("numtrees", "mtry", "OOB_mse")
# library
library(scatterplot3d)



scatterplot3d(x = data$numtrees,
              y = data$mtry,
              z = data$OOB_mse,
              xlab = "Number of Trees",
              ylab = "Number of Features",
              zlab = "OOB MSE",
              main = "Number of Trees & Number of Features vs. OOB MSE",
              angle = 60,
              pch = 16,
              color = ifelse((data$OOB_mse < 769.72), "red", "#316D9E"),
              type = "h")

sorted <- sort(data$OOB_mse)
sorted

