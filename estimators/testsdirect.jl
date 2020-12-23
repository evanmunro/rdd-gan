using CSV, DataFrames
includet("directrdd.jl")

data = DataFrame(CSV.File("../data/cleaned/lee.csv"))
data = first(data, 500)

localData(data.y, data.x)
