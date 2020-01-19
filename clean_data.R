
#Lee Data
data <- read.csv("~/Documents/Github/rdd-gan/data/lee/davidlee.csv") 
write.csv(data.frame(x=data$x,y=data$y1),file="data/cleaned/lee.csv",row.names=F) 

