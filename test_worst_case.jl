using CSV, DataFrames
include("estimators.jl")
include("worst_case_fn.jl")
include("rdddata.jl") 

M = 14.28
data = first(DataFrame(CSV.File("data/cleaned/lee.csv")), 500)
data.x = round.(data.x, digits=3)
data = RDData(data.y, data.x)

tau_opt, weights_opt = optrdd(data, data.sig2, 14.28, true)
println(tau_opt)
tau_llr, weights_llr = llrrdd(data, triangular_kernel, compute_opt_bw(data), true)

γ0_opt, γ1_opt = get_gamma(data, weights_opt)
plot([data.xx0..., data.xx1...], [γ0_opt..., γ1_opt...])
γ0_llr, γ1_llr = get_gamma(data, weights_llr)
#check weights
plot(data.xx0, [γ0_opt γ0_llr],  label=nothing)

#calculate worst case fn

result = worst_case_mu(data, γ0_opt, γ1_opt, 14.28)
result = worst_case_mu(data, γ0_llr, γ1_llr, 14.28)
#plot worst case fn
plot([data.xx0...,  data.xx1...], result.max_mu)
plot!([data.xx0...,  data.xx1...], [data.mu0..., data.mu1...], seriestype=:scatter)
