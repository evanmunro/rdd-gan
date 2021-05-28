using RCall, LinearAlgebra, GLM
using Suppressor
include("rdddata.jl")

R"library(optrdd)"

function estimate_B2(x, y)
	@rput x; @rput y;
	@suppress R"""
		B = RDHonest::NPR_MROT.fit(RDHonest::RDData(data.frame(y=y, x=x), cutoff=0))
	"""
	@rget B
end
function compute_opt_bw_cv(x ,y)
	c=0
	@rput x; @rput y; @rput c
	R"""
		bw_selection = rdrobust::rdbwselect_2014(y, x, c=c, bwselect="CV")
		opt_bw = bw_selection[["bws"]]
	"""
	@rget opt_bw
	return opt_bw[1]
end

function optrdd_R(x, y, ret_weights=false)
	c = 0; w = x .>= 0
	B = estimate_B2(x, y)
	sigma2 = estimate_σ2(x, y)
	@rput x; @rput y; @rput B; @rput c; @rput w; @rput sigma2;
	@suppress R"""
		out = optrdd(X=x, W=w, Y=y, estimation.point=c,
					max.second.derivative=B, sigma.sq=sigma2,
					optimizer="mosek", use.spline=FALSE, try.elnet.for.sigma.sq=FALSE,
					verbose=FALSE)
		tau = out[["tau.hat"]]
	"""
	@rget tau

	if ret_weights
		R"""gamma = out[["gamma"]]"""
		@rget gamma
		return tau, gamma
	end

	return tau
end

triangular_kernel(u) = (abs(u) <= 1)*(1-abs(u))

# NOTE: only implemented for triangular kernel right now
function compute_opt_bw(x, y)
	c=0
	@rput x; @rput y; @rput c
	R"""
		bw_selection = rdrobust::rdbwselect_2014(y, x, c=c, bwselect="IK")
		opt_bw = bw_selection[["bws"]]
	"""
	@rget opt_bw
	return opt_bw[1]
end


function llrrdd(x, y, kernel, bw, ret_weights=false)
	c=0
	n = length(x)
	@rput y; @rput x; @rput c; @rput bw
	R"""
		model <- rdd::RDestimate(y~x, data.frame(y=y, x=x), cutpoint=c, bw=bw)
		est <- model$est[1]
		se <-
	"""
	@rget est
	if ret_weights
		Δ = x .- c; w = x .>= c
		X = [ones(n) w min.(Δ, 0) max.(Δ, 0)]
		K = Diagonal(kernel.(abs.(Δ)./bw))
		H = inv(X'*K*X)*X'*K
		weights = H[2, :]
		return est, weights
	end
	return est
end
