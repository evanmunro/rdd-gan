using CSV, DataFrames
includet("directrdd.jl")

data = DataFrame(CSV.File("../data/cleaned/lee.csv"))
data = first(data, 500)

#localData(data.y, data.x)
y  = data.y ; x= data.x
y1 = y[x.>0]
x1 = x[x.>0]

#dWeightsConstant(y1, x1)

rddConstant(data.y, data.x)

rddLasso(data.y, data.x)
rddPoly(data.y, data.x)
