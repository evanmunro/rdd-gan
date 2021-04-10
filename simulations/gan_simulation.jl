using Feather, Statistics, StatsBase, CSV, DataFrames
includet("../../MinMaxRDD/estimators.jl")
includet("../../MinMaxRDD/optrdd.jl")
includet("../estimators/directrdd.jl")
#include("adaptiveminmax.jl")

dataset = "lee"
digits = nothing
gen_path = string("../data/generated/", dataset, "_generated.feather")
real_path = string("../data/cleaned/", dataset, ".csv")



df_real = DataFrame!(CSV.File(real_path))
n1  = sum(df_real.x .> 0)
n0 = sum(df_real.x .<= 0)
df = Feather.read(gen_path)
df.x = convert(Array{Float64}, df.x)
df.y = convert(Array{Float64}, df.y)
t_gt = llrrdd(df.x, df.y, triangular_kernel, compute_opt_bw(df.x, df.y), false)
println("ground truth: ", t_gt)

if !(digits==nothing)
    df.x = round.(df.x, digits=digits)
end

df = df[df.x.!=0, :]
idx1   = collect(1:nrow(df))[df.x .>0]
idx0 = collect(1:nrow(df))[df.x .<=0]

sample_data(df) = vcat(df[sample(idx1, n1, replace=false), :],
                       df[sample(idx0, n0, replace=false), :] )

#println("tau: ", result.τ)
estimators = ["LLR", "CV-LM", "CV-NEW"]
nestimators = length(estimators)
nsims = 100
t_est = zeros(nsims, nestimators)
#covg_results = zeros(nsims, nestimators)

for i in 1:nsims
    data = sample_data(df)
    t_est[i, 1] = llrrdd(data.x, data.y, triangular_kernel, compute_opt_bw(data.x, data.y), false)
    #t_est[i, 1] = rddCV(data.y, data.x)
    t_est[i, 2] = llrrdd(data.x, data.y, triangular_kernel, compute_opt_bw_cv(data.x, data.y), false)
    t_est[i, 3] = rddCVSingle(data.y, data.x)
    println(t_est[i,:])
end

valid = t_est[:, 3] .!= 0.0
#t_trim = t_est[t_est[:, 3] .< maximum(t_est[:, 2]), :]
println("llr RMSE: ", sqrt(mean((t_est[valid, 1].- t_gt).^2)))
println("optrdd RMSE: ",sqrt(mean((t_est[valid, 2] .- t_gt).^2)))
println("adaptrdd RMSE: ",sqrt(mean((t_est[valid, 3] .- t_gt).^2)))

println("llr bias: ", mean(t_est[valid, 1].- t_gt))
println("optrdd bias: ", mean(t_est[valid, 2] .- t_gt))
println("adaptrdd bias: ", mean(t_est[valid, 3] .- t_gt))

println("llr SD:  ", std(t_est[valid, 1]))
println("optrdd SD: ", std(t_est[valid, 2]))
println("adaptrdd SD: ", std(t_est[valid, 3]))

#println("adaptrdd RMSE TRIM: ",sqrt(mean((t_trim[:, 3] .- t_gt).^2)))
