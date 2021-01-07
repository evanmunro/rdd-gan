using CSV, DataFrames
includet("directrdd.jl")

data = DataFrame(CSV.File("../data/cleaned/mats_math.csv"))
data = first(data, 500)

rddCV(data.y, data.x)

#localData(data.y, data.x)
y  = data.y ; x= data.x
y1 = y[x.>0]
x1 = x[x.>0]
y0 = y[x.<0]
x0 = x[x.<0]
optimalWeights(y1, x1, 0.15, 1)

optimalWeights(y0, -1 .* x0, 0.15, 1)
#dWeightsConstant(y1, x1)

rddConstant(data.y, data.x)

rddLasso(data.y, data.x)
rddPoly(data.y, data.x)
