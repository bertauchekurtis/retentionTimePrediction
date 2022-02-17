# make a dataset for AUTO RT
# kurtis bertauche

library(stringr)

data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")

data$Peptide.Sequence2 <- str_replace_all(data$Peptide.Sequence2, "m", "1")
data$Peptide.Sequence2 <- str_replace_all(data$Peptide.Sequence2, "s", "2")
data$Peptide.Sequence2 <- str_replace_all(data$Peptide.Sequence2, "t", "3")
data$Peptide.Sequence2 <- str_replace_all(data$Peptide.Sequence2, "y", "4")

names(data) <- c("x", "y")

set.seed(37)
setAssignments <- sample(1:2, size = nrow(data), prob = c(0.8, 0.2), replace = TRUE)
trainingData <- data[setAssignments == 1,]
testingData <- data[setAssignments == 2,]

write.table(trainingData, 
            file = "C:/Users/Kurtis/Downloads/trainingSet.tsv", 
            row.names=FALSE, 
            quote = FALSE, 
            sep="\t")

write.table(testingData, 
            file = "C:/Users/Kurtis/Downloads/testing.tsv", 
            row.names=FALSE, 
            quote = FALSE, 
            sep="\t")

