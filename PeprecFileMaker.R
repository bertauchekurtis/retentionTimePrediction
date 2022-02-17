# kurtis bertauche
# peprec file maker

library(stringr)

data <- read.csv(file = "C:/Users/Kurtis/Desktop/Research/data/RetentionTime_HCD_Marx2013_SuppT3.csv")

# function to parse and create modification list
parser = function(string)
{
  resultingString = ""
  
  for(position in 1:nchar(string))
  {
    if(substr(string, position, position) == "m")
    {
      resultingString <- paste(resultingString, paste(position, "Oxidation", sep= "|"), sep = "|")
    }
    else if ((substr(string, position, position) == "s")|(substr(string, position, position) == "t")|(substr(string, position, position) == "y"))
    {
      resultingString <- paste(resultingString, paste(position, "Phospho ", sep= "|"), sep = "|") 
    }
  }
  sub('.','',resultingString)
}

# apply function
data$modifications <- sapply(data$Peptide.Sequence2, parser)

# rename columns
colnames(data) <- c("seq", "tr", "modifications")

# reorer columns
data <- data[, c(1, 3, 2)]

# change sequence data to match file format
data$seq <- str_replace_all(data$seq, "m", "M")
data$seq <- str_replace_all(data$seq, "s", "S")
data$seq <- str_replace_all(data$seq, "t", "T")
data$seq <- str_replace_all(data$seq, "y", "Y")

# split data
set.seed(37)
setAssignments <- sample(1:2, size = nrow(data), prob = c(0.8, 0.2), replace = TRUE)
trainingData <- data[setAssignments == 1,]
testingData <- data[setAssignments == 2,]

# export data
write.table(trainingData, 
            file = "C:/Users/Kurtis/Downloads/PEPREC_trainingSet.csv", 
            row.names=FALSE, 
            quote = FALSE, 
            sep=",")

write.table(testingData, 
            file = "C:/Users/Kurtis/Downloads/PEPREC_testing.csv", 
            row.names=FALSE, 
            quote = FALSE, 
            sep=",")
