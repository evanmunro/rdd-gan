
library(feather)

x = runif(100000, min=-1, max=1)
epsilon = rnorm(100000)
ysmooth = (x<0)*(2*x^3) + (x>=0)*(2*x^3 + 0.5 - 0.5*x^2) +epsilon
ycurved = (x<0)*(10*((x+0.5)^3 - 0.5^3) + 0.2*((x-0.7)^2 - (-0.7)^2) + 0.1*((x+1)^2-1^2)) +
          (x>=0)*(0.5+(40*((x-0.2)^3-(-0.2)^3)-80*((x-0.2)^4-(-0.2)^4)+38*((x-0.2)^5-(-0.2)^5))) +
          epsilon

data_smooth=data.frame(x=x, y=ysmooth)

data_curved=data.frame(x=x,y=ycurved)
write_feather(data_smooth, "data/generated/mc_smooth.feather")
write_feather(data_curved, "data/generated/mc_curved.feather")
