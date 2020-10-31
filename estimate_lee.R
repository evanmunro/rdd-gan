data <- read.csv("data/cleaned/lee.csv")
data <- data[abs(data$x) < 0.18,]
X = data$x
Y= data$y
c= 0 
M= RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))
sigma = 
W = as.numeric(X>c)
model <-  optrdd::optrdd(X,Y,W,
                         max.second.derivative=M,verbose=F,estimation.point=c,
                         sigma.sq = 0.0118858279376246, 
                         try.elnet.for.sigma.sq=F, use.spline=F, optimizer='mosek')
