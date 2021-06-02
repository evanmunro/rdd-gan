using Plots, DataFrames, CSV
includet("simulations/ikhonest.jl")

data = DataFrame!(CSV.File("data/cleaned/lee.csv"))


model = RDDIK(data.x, data.y)
bw = model.fitted_bandwidth
ik_weights = calculate_weights(data.y, data.x, bw)

xR = RunningVariable(data.x; cutoff= 0.0)
rd = RDData(data.y, xR)
model = fit(ImbensWagerOptRD(B=14.28, solver=Mosek.Optimizer), rd.ZsR, rd.Ys)
iw_weights = model.weights

sorted_data_perm = sortperm(data.x)

weight_data = hcat(data.x[sorted_data_perm], ik_weights[sorted_data_perm], iw_weights[sorted_data_perm])
unique_data = unique(weight_data, dims=1)


plot(unique_data[:, 1], unique_data[:, 2:3], xrange= [ -0.3, 0.3], labels= ["IK" "IW"],
            xlabel="Running Variable", ylabel="Weights")

savefig("lee_weights.pdf")
