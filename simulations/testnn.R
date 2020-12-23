source("estimators.R")


data <- read.csv("../data/cleaned/lee.csv")
X = data$x
Y= data$y

rddNN(X,Y, 0)