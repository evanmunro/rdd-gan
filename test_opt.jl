using CSV, DataFrames, Plots
include("optrdd.jl")
include("estimators.jl")

data = DataFrame(CSV.File("data/cleaned/lee.csv"))
#data = data[abs.(data.x).<1 , :]
data = data[abs.(data.x).<0.18 , :]
#data = data[abs.(data.x).<50 , :]
@time test = optrdd(data.x, data.y)
println("Julia Estimate: ", test.τ)
@time tau_opt, weights_opt = optrdd_R(data.x, data.y, true)
println("R Estimate: ", tau_opt)
plot(data.x, [test.γ, weights_opt], seriestype=:line, label=["Julia Weights" "R Weights"]) #xlim=[-0.15, 0.15])
