data <- read.csv("data/cleaned/lee.csv")
X = data$x
Y = data$y 
M=RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))

fit = rdd::RDestimate(y~x,data,model=TRUE,frame=TRUE)
bw= fit$bw[[1]]
wts.ik = rdd::kernelwts(data$x,0,bw)

W = as.numeric(data$x>0)
model <-  optrdd::optrdd(data$x,data$y,W,
                         max.second.derivative=M,verbose=F,estimation.point=0,
                         try.elnet.for.sigma.sq=T,optimizer='mosek')
wts.iw <- abs(model$gamma)/sum(abs(model$gamma))

model <- RDHonest::RDHonest(y~x,data=gen,cutoff=0,M=14.28,opt.criterion="MSE")