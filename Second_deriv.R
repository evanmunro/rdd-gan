r1 <- unname(stats::lm(data$y ~ 0 + outer(data$x, 0:4, "^"))$coefficients)
f2 <- function(x) abs(2*r1[3]+6*x*r1[4]+12*x^2*r1[5])

d2 = f2(data$x)

global_M <- function(df) { 
  model <- stats::lm(df$y ~ 0 + outer(df$x, 0:4, "^"))
  pdf()
  plot(df$x, predict(model), main="poly")
  Sys.sleep(4)
  r1 <- unname(model$coefficients)
  f2 <- function(x) abs(2*r1[3]+6*x*r1[4]+12*x^2*r1[5])
  plot(df$x, f2(df$x), main="d2")
  return(f2(df$x))
}

check_M <- function(path) { 
  data <- read.csv(path) 
  da <- data[data$x>0, ]
  db <- data[data$x<=0, ]
  
  rd <- rddtools::rdd_data(x=data$x, y=data$y, cutpoint=0)
  bw <- rddtools::rdd_bw_ik(rd)
  print("global")
  a2 = global_M(da) 
  print(max(a2))
  b2 = global_M(db) 
  print(max(b2))
  
  print("ik")
  print(max(a2[abs(da$x)<bw]))
  print(max(b2[abs(db$x)<bw]))
}
