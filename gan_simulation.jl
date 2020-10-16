using Feather, Statistics, StatsBase
include("estimators.jl")
include("optrdd.jl")

datapath = "data/generated/lee_generated.feather"
df = Feather.read(datapath)
#ya = df.y[0 .<= df.x .<= 0.00001]
#yb = df.y[-0.00001 .<= df.x .<= 0]
#t_gt = mean(ya) - mean(yb)
n = 6558
B = 14.28108

t_gt = llrrdd(df.x, df.y, triangular_kernel, compute_opt_bw(df.x, df.y), false)

nsims = 200
t_est = zeros(200,2)
idx1   = collect(1:nrow(df))[df.x .>0]
idx0 = collect(1:nrow(df))[df.x .<=0]
for i in 1:nsims
    data = vcat(df[sample(idx1, 3818, replace=false), :],
             df[sample(idx0, 2740, replace=false), :] )
    #data = OptRDData(data.y, data.x)
    t_est[i, 1] = llrrdd(data.x, data.y, triangular_kernel, compute_opt_bw(data.x, data.y), false)
    t_est[i, 2] = optrdd(data.x, data.y, 0.01915)
end

println("llr RMSE: ", sqrt(mean((t_est[:, 1].- t_gt).^2)))
println("optrdd RMSE: ",sqrt(mean((t_est[:, 2] .- t_gt).^2)))

println("llr bias: ", mean(t_est[:, 1].- t_gt))
println("optrdd bias: ", mean(t_est[:, 2] .- t_gt))

println("llr SD:  ", std(t_est[:, 1]))
println("optrdd SD: ", std(t_est[:, 2]))
