
#Lee Data
data <- read.csv("~/Documents/Github/rdd-gan/data/lee/davidlee.csv") 
write.csv(data.frame(x=data$x,y=data$y1),file="data/cleaned/lee.csv",row.names=F) 

#Matsudaira data 
library(R.matlab) 
raw.path = "~/Dropbox/rdd_gan/data/matsudaira/matsudaira_11apr2.mat"
data <- readMat(raw.path) 

math.data <- data.frame(y=data$zmscr02,x=data$mdcut01,t=data$ssatyn01)
read.data <- data.frame(y=data$zrscr02,x=data$rdcut01,t=data$ssatyn01)

write.csv(math.data,file="~/Documents/Github/rdd-gan/data/cleaned/m_math.csv",row.names=F)
write.csv(read.data,file="~/Documents/Github/rdd-gan/data/cleaned/m_read.csv",row.names=F)
#outcome is either zmscr02 or zrscr02 which is 2002 math and reading score 
#summer school attendance is ssatyn01 (that is treatment) 

#2000 math and reading score is mmscr00 and mrscr00 

#pass fail cutoff for math and reading is mdcut01 and rdcut01 

#All observations are grade 5 

#Jacob-Lefgren data 
library(readstata13)
raw.data <- read.dta13("~/Dropbox/rdd_gan/data/jacob_lefgren/imbens_small.dta")
math.data <- data.frame(y=raw.data$rmath1,z1=raw.data$mtge0-2.8,z2 = raw.data$rdge0-2.8)
math.data$x <- pmin(math.data$z1,math.data$z2)
math.data$w <- math.data$x<0 
write.csv(math.data,file="~/Documents/Github/rdd-gan/data/cleaned/jl_math.csv",row.names=F)

