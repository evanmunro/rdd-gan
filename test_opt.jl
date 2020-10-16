using CSV, DataFrames, Plots

include("rdddata.jl")
include("optrdd.jl")
data = DataFrame(CSV.File("data/cleaned/lee.csv"))
data = data[abs.(data.x).<0.18, :]
@time test = optrdd_jl(data.x, data.y)
println("Julia Estimate: ", test.τ)
include("estimators.jl")
@time tau_opt, weights_opt = optrdd(data.x, data.y, true)
println("R Estimate: ", tau_opt)
plot(data.x, [test.γ, weights_opt], seriestype=:line, label=["Julia Weights" "R Weights"], xlim=[-0.15, 0.15])
