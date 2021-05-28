using CSV, DataFrames
includet("simulations/ikhonest.jl")
data = DataFrame!(CSV.File("data/cleaned/lee.csv"))
