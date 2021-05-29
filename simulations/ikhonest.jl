using JuMP, GLM, RegressionDiscontinuity, MosekTools, LinearAlgebra
using Roots: fzero
using Loess: loess, predict
using FiniteDifferences
using DataFrames, Suppressor

struct DiscretizedRDData
	y::Array{Float64, 1}
	x::Array{Float64, 1}
	xx::Array{Float64, 1}
	weights::Array{Int, 1}
	ixc::Int
	h::Array{Float64, 1}
	xmap::Array{Int, 1}
    grid::Array{Float64, 1}
	σ²::Float64
end

function smooth_mean(x ,y)
	model = loess(x, y)
	return f(a) = predict(model, a)
end

function fit_bounds_l(x, y)

	xcut = quantile(abs.(x), 0.8)
	println("xcut::", xcut)
	y = y[abs.(x) .< xcut]; x = x[abs.(x) .< xcut]

	μ = smooth_mean(x,y)
	us = range(extrema(x)...; step = 0.1)
	f1 = zeros(length(us))
	f2 = zeros(length(us))
	for i in 2:(length(us) -1)
		f1[i] = central_fdm(5, 1)(μ, us[i])
		f2[i] = central_fdm(5, 2)(μ, us[i])
	end
	return (M1 = maximum(f1), M2 = maximum(f2))
end

function estimate_M(data::DiscretizedRDData)
	w = data.x .> 0
	y1 = data.y[w.==1]; x1 = data.x[w.==1]; y0 = data.y[w.==0]; x0 = data.x[w.==0]
	boundsa = fit_bounds_l(x1, y1)
	boundsb = fit_bounds_l(x0, y0)
	return (M1 = max(boundsa.M1, boundsb.M1), M2=max(boundsa.M2, boundsb.M2))
end

function fit_bounds(x, y)
	xcut = quantile(abs.(x), 1.0)
	println("xcut::", xcut)
	y = y[abs.(x) .< xcut]; x = x[abs.(x) .< xcut]
	model = lm(@formula(y ~ x + x^2 + x^3 + x^4), DataFrame(y=y, x=x))
	beta = coef(model)
	f1(x) = beta[2] .+ 2*beta[3].*x + 3*beta[4].*x.^2 + 4*beta[5].*x.^3
    f2(x) = 2*beta[3] .+ 6*beta[4].*x .+ 12*beta[5]*x.^2
    f2_est = abs.(f2(x))
	f1_est = abs.(f1(x))
	println("f1: ", maximum(f1_est))
	println("f2: ", maximum(f2_est))
    return (M1 = maximum(f1_est), M2 = maximum(f2_est))
end

function add_constraint!(model, f, h, M1, M2)
    d  = length(h) + 1
    for i in 3:d
        @constraint(model, -M2 <= (f[i] - 2*f[i-1] + f[i-2])/(h[i-1]*h[i-2])<= M2)
    end
	for i in 2:d
         @constraint(model, 0 <= (f[i] - f[i-1])/h[i-1] <= M1)
    end
end

function DiscretizedRDData(y, x, num_buckets)
	xmax = maximum(x)
 	xmin = minimum(x)
	h = (xmax - xmin)/num_buckets
	grid = (xmin-h/2):h:(xmax+h)
	grid = vcat(grid, repeat([0], 2 - sum(grid.==0)))
	sort!(grid)
	xx = midpoints(grid)
	h  = xx[2:length(xx)] .- xx[1:(length(xx)-1)]
	ixc = argmin(abs.(xx))
	binned = fit(Histogram, x, grid)
	weights = binned.weights
	xmap = StatsBase.binindex.(Ref(binned), x)
	sigma2 = estimate_σ2(x, y)
	return DiscretizedRDData(y, x, midpoints(grid), weights, ixc, h, xmap, grid, sigma2)
end

function estimate_σ2(x, y)
	data = DataFrame(y=y, x=x, w= x.>0)
    yhat = GLM.predict(lm(@formula(y~x*w),data))
    σ2 = mean((data.y - yhat).^2) * length(data.x) / (length(data.x) - 4)
	return σ2
end

function test_worst_case(x, y)
	x = RunningVariable(x; cutoff= 0.0)
	rd = RDData(y, x)
	model = fit(ImbensWagerOptRD(B=14.28, solver=Mosek.Optimizer), rd.ZsR, rd.Ys)
	println("tau_est", model.tau_est)
	println("se: ", model.se_est)
	println("bias: ", worst_case_mu(x, y, model.weights).max_abs_bias)
end

function worst_bias_symmetric(x, y, γ)
	model =  Model(Mosek.Optimizer)
	println("sum gamma: ", sum(γ))

	data  = DiscretizedRDData(y, x, 2000)
	bounds = estimate_M(x,y)
	M1 = bounds.M1
	M2 = bounds.M2 #estimate_M(x,y)
	println(M1)
	X = data.xx; d=length(X); h = data.h; n= data.weights; ixc = data.ixc
	j = data.xmap
	@variable(model, mu[1:d])
	add_constraint!(model, mu, h, M1, M2)
	@constraint(model, mu[ixc] == 0 )
    @constraint(model, (mu[ixc+1] - mu[ixc])/h[ixc] == 0)
	@objective(model, Max, (sum(γ[i]*mu[j[i]] for i in 1:sum(n))))
	optimize!(model)
    max_bias = objective_value(model)
	return max_bias

end
function worst_case_mu(x, y, γ, γ0 = 0.0)
    model =  Model(Mosek.Optimizer)
    data  = DiscretizedRDData(y, x, 2000)
	println("sigma: ", data.σ²)
	bounds = estimate_M(data)
	M1 = bounds.M1
	M2 = bounds.M2
	println("M1 bound", M1)
	println("M2 Bound", M2)
    #M2 = 14.28#estimate_M(x,y)
	#M=1200
    X = data.xx; d=length(X); h = data.h; n= data.weights; ixc = data.ixc
    @variable(model, mu0[1:d])
    @variable(model, mu1[1:d])
    c = 0.0; W = X .> c; j = data.xmap
	#@constraint(model, mu0[ixc]==0)
	#@constraint(model, mu1[ixc]==0)
	#@constraint(model, (mu0[ixc+1] - mu0[ixc])/h[ixc] == 0)
   	#@constraint(model, (mu1[ixc+1] - mu1[ixc])/h[ixc] == 0)

	add_constraint!(model, mu0, h, M1, M2)
	add_constraint!(model, mu1, h, M1, M2)

    # maximize the absolute bias
    @objective(model, Max, (sum(γ[i]*mu1[j[i]]*W[j[i]] for i in 1:sum(n)) +
                            sum(γ[i]*mu0[j[i]]*(1-W[j[i]]) for i in 1:sum(n))
                            + γ0 - mu1[ixc] + mu0[ixc]))

    optimize!(model)
    max_bias = objective_value(model)

    max_mu = vcat(value.(mu0)[W.==0], value.(mu1)[W.==1])

    @objective(model, Min, (sum(γ[i]*mu1[j[i]]*W[j[i]] for i in 1:sum(n)) +
                            sum(γ[i]*mu0[j[i]]*(1-W[j[i]]) for i in 1:sum(n))
                            + γ0 - mu1[ixc] + mu0[ixc]))

    optimize!(model)
    min_bias = objective_value(model)
    min_mu = vcat(value.(mu0)[W.==0], value.(mu1)[W.==1])

    max_abs_bias = max(abs(max_bias), abs(min_bias))

    return (max_abs_bias = max_abs_bias,
            max_bias = max_bias,
            min_bias = min_bias,
            max_mu = max_mu,
            min_mu = min_mu,
            xx = X)
end

function bias_adjusted_gaussian_ci(se; maxbias=0.0, level=0.95)
    rel_bias = maxbias/se
    zz = fzero( z-> cdf(Normal(), rel_bias-z) + cdf(Normal(), -rel_bias-z) +  level -1,
        0, rel_bias - quantile(Normal(),(1- level)/2.1))
    zz*se
end

function RDDIK(x, y)
	x = RunningVariable(x; cutoff= 0.0)
	data = RegressionDiscontinuity.RDData(y, x)
	model = fit(NaiveLocalLinearRD(kernel = SymTriangularDist(),
		bandwidth = ImbensKalyanaraman()), data.ZsR, data.Ys)
	return model
end

triangular_kernel(u) = (abs(u) <= 1)*(1-abs(u))

function calculate_weights(y, x, bw)
	n = length(x)
	Δ = x; w = x .>= 0.0
	X = [ones(n) w min.(Δ, 0) max.(Δ, 0)]
	K = Diagonal(triangular_kernel.(abs.(Δ)./bw))
	H = inv(X'*K*X)*X'*K
	weights = H[2, :]
end

function IKHonest(y, x)
	model = RDDIK(x, y)
	bw = model.fitted_bandwidth
	se_est = model.se_est
	weights = calculate_weights(y, x, bw)
	max_bias = 0.0
	@suppress_out begin
		max_bias = worst_case_mu(x, y, weights).max_abs_bias
	end
	#println("bias: ", max_bias)
	ci = bias_adjusted_gaussian_ci(se_est; maxbias = max_bias)
	return (ate = model.tau_est,
			se = se_est,
			lower = model.tau_est - ci,
			upper = model.tau_est + ci
			)
end
