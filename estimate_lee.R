data <- head(read.csv("data/cleaned/lee.csv"), 500)

X = data$x
Y= data$y
c= 0 
M= RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))
W = as.numeric(X>c)
model <-  optrdd::optrdd(X,Y,W,
                         max.second.derivative=M,verbose=F,estimation.point=c,
                         try.elnet.for.sigma.sq=T,optimizer='mosek')
