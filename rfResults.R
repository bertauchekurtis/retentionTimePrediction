# kurtis bertauche
# rf visualization

# load data
data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/resultsRFFixed.csv")

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
              color = ifelse((data$OOB_mse < 26.41394), "red", "#316D9E"),
              type = "h")

sorted <- sort(data$OOB_mse)
sorted

