
global_M <- function(df, name) { 
  fname = paste0(name,".pdf")
  model <- stats::lm(df$y ~ 0 + outer(df$x, 0:4, "^"))
  pdf(file=paste0(name,"smooth.pdf"), width=4, height=4)
  plot(df$x, predict(model), main="Poly Regression",xlab="x", ylab="Smoothed")
  dev.off()
  r1 <- unname(model$coefficients)
  f2 <- function(x) abs(2*r1[3]+6*x*r1[4]+12*x^2*r1[5])
  f1 <- function(x) abs(4*x^3*r1[5] + 3*x^2*r1[4] + 2*x*r1[3] + r1[2]) 
  pdf(file=paste0(name,"d2.pdf"), width=4, height=4)
  plot(df$x, f2(df$x), main="Second Derivative", xlab="x", ylab="D2")
  dev.off() 
  return(list(f1(df$x), f2(df$x)))
}

check_M <- function(name) { 
  data <- read.csv(paste0("data/cleaned/",name,".csv") )
  da <- data[data$x>0, ]
  db <- data[data$x<=0, ]
  
  rd <- rddtools::rdd_data(x=data$x, y=data$y, cutpoint=0)
  bw <- rddtools::rdd_bw_ik(rd)
  print("global")
  result = global_M(da, paste0(name,"a")) 
  a1 = result[[1]]
  a2 = result[[2]]
  print(max(a1))
  print(max(a2))
  result = global_M(db, paste0(name,"b"))
  b1 = result[[1]]
  b2 = result[[2]]
  print(max(b1))
  print(max(b2))
  
  print("ik")
  print(max(a2[abs(da$x)<bw]))
  print(max(b2[abs(db$x)<bw]))
}

check_M("lee")
check_M("jl_math")
check_M("m_math")