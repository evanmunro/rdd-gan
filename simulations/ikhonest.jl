using JuMP, GLM, RegressionDiscontinuity, MosekTools, LinearAlgebra
using Roots: fzero

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

function estimate_M(x, y)
	model = lm(@formula(y ~ x + x^2 + x^3 + x^4), DataFrame(y=y, x=x))
	beta = coef(model)
    f2(x) = 2*beta[3] .+ 6*beta[4].*x .+ 12*beta[5]*x.^2
    f2_est = abs.(f2(x))
    return maximum(f2_est)
end

function add_constraint!(model, f, h, M2)
    d  = length(h) + 1
    for i in 3:d
        @constraint(model, M2 <= (f[i] - 2*f[i-1] + f[i-2])/(h[i-1]*h[i-2])<= M2)
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

function worst_case_mu(x, y, γ, γ0 = 0.0)
    model =  Model(Mosek.Optimizer)
	println(typeof(y))
	println(typeof(x))
    data  = DiscretizedRDData(y, x, 2000)
    M = estimate_M(x,y)
	#M=1200
    X = data.xx; d=length(X); h = data.h; n= data.weights; ixc = data.ixc
    @variable(model, mu0[1:d])
    @variable(model, mu1[1:d])
    c = 0.0; W = X .> c; j = data.xmap

    add_constraint!(model, mu0, h, M)
    add_constraint!(model, mu1, h, M)
    ## add these constraints for model to have a solution
    @constraint(model, mu0[ixc] == 0 )
    @constraint(model, mu1[ixc] == 0 )
    @constraint(model, (mu0[ixc+1] - mu0[ixc])/h[ixc] == 0)
    @constraint(model, (mu1[ixc+1] - mu1[ixc])/h[ixc] == 0)

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
	data = RDData(y, x)
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
	max_bias = worst_case_mu(x, y, weights).max_abs_bias
	println(max_bias)
	ci = bias_adjusted_gaussian_ci(se_est; maxbias = max_bias)
	return (ate = model.tau_est,
			se = se_est,
			lower = model.tau_est - ci,
			upper = model.tau_est + ci
			)
end
