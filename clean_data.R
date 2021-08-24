
#Lee Data
data <- read.csv("~/Documents/Github/rdd-gan/data/lee/davidlee.csv")
write.csv(data.frame(x=data$x,y=data$y1),file="data/cleaned/lee.csv",row.names=F)

#Matsudaira data
library(R.matlab)
raw.path = "~/Dropbox/rdd_gan/data/matsudaira/matsudaira_11apr2.mat"
data <- readMat(raw.path)

math.data <- data.frame(y=data$zmscr02, x=pmin(data$mdcut01, data$rdcut01), t=data$ssatyn01)
read.data <- data.frame(y=data$zrscr02, x=pmin(data$mdcut01, data$rdcut01), t=data$ssatyn01)

write.csv(math.data,file="~/Documents/Github/rdd-gan/data/cleaned/mats_math.csv",row.names=F)
write.csv(read.data,file="~/Documents/Github/rdd-gan/data/cleaned/mats_read.csv",row.names=F)
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


#Meyersson (2014) data on Turkey elections
#x is min -100, max 100 , y is [0, 100]
data <- read.csv("data/meyersson_raw.csv") 
data <- data.frame(y = data$Y, x= data$X)
write.csv(data, "data/cleaned/meyersson.csv")
#Senate Incumbency data on US elections
#x is min -100, max 100 , y is [0, 100]
data <- read.csv("data/senate_raw.csv")
data <- data.frame(y = data$demvoteshfor2, x= data$demmv)
data <- data[!is.na(data$y), ]
write.csv(data, "data/cleaned/senate.csv")
#Incumbency curse on Brazil electoins
#x is min -100, max 100 , y is [-100, 100]
data <- read.csv("data/brazil_raw.csv")
data <- data.frame(y = data$mv_incpartyfor1, x= data$mv_incparty)
data <- data[!is.na(data$y), ]
data <- data[!is.na(data$x), ]
write.csv(data, "data/cleaned/brazil.csv")
