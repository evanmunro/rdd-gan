#0.05 to 0.7
using CSV, DataFrames, Plots, Feather
includet("directrdd.jl")
includet("../../MinMaxRDD/estimators.jl")

dataset = "lee"
digits = nothing
gen_path = string("../data/generated/", dataset, "_generated.feather")

data = DataFrame(CSV.File("../data/cleaned/lee.csv"))

n1 = 3818
n0 = 2740
df = Feather.read(gen_path)
df.x = convert(Array{Float64}, df.x)
df.y = convert(Array{Float64}, df.y)
df = df[df.x.!=0, :]
idx1   = collect(1:nrow(df))[df.x .>0]
idx0 = collect(1:nrow(df))[df.x .<=0]

sample_data(df) = vcat(df[sample(idx1, n1, replace=false), :],
                       df[sample(idx0, n0, replace=false), :] )

γs = 0.05:0.05:0.9
cvCriterionlee = zeros(length(γs))
cvGan = zeros(length(γs))
sdata = sample_data(df)
for i in 1:length(γs)
    cvCriterionlee[i] = rddCVObjective(γs[i], data.y, data.x)
    cvGan[i] = rddCVObjective(γs[i], sdata.y, sdata.x)
end

cvGan = cvGan.*mean(cvCriterionlee)./mean(cvGan)
plot(γs, [cvCriterionlee cvGan], labels=["Lee Data" "GAN Sample"])
savefig("CV_criterion.pdf")

t_gt = 0.09355633346702093

nsims = 100
bw = zeros(nsims)
bwcv = zeros(nsims)

mse = zeros(length(γs))
est = zeros(length(γs))
bias² = zeros(length(γs))
var = zeros(length(γs))
for j in 1:nsims
    data = sample_data(df)
    bw[j] = compute_opt_bw(data.x, data.y)
    bwcv[j] = cvBWSelect(data.y, data.x)
    println(bwcv[j])
    println(bw[j])
    for i in 1:length(γs)
        γ = γs[i]
        τ =rddCV(data.y, data.x, γ)
        est[i] += τ/nsims
        mse[i] += (τ - t_gt)^2/nsims
    end
end

rmse = sqrt.(mse)
bias² = (est .- t_gt).^2
var = mse .- bias²

var = var./mean(var).*mean(rmse)
bias² = bias²./mean(bias²).*mean(rmse)

plot(γs, [rmse bias² var], labels=["RMSE" "Bias" "Var" "Mean_IK_BW"])
vline!([mean(bw) mean(bwcv)], labels=["Mean_IK_BW" "Mean_CV_BW"])
savefig("rmseplot.pdf")


histogram(bw, xlim=[0.1,1], label="IK_BW")
vline!([γs[argmin(rmse)]], label="Optimal")
savefig("bw_IK.pdf")
histogram(bwcv, xlim=[0.0, 1], label="CV_BW")
vline!([γs[argmin(rmse)]], label="Optimal")
savefig("bw_CV.pdf")

plot(bw, bwcv, seriestype=:scatter, xlabel="IK_BW", ylabel="CV_BW")
savefig("bw_scatter.pdf")
