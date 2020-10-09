using RCall, LinearAlgebra
using Suppressor
include("rdddata.jl")

R"library(optrdd)"

function optrdd(data::RDData, sigma_sq, B, ret_weights=false)
	x = data.x; y = data.y; c = 0; w = x .>= 0
	@rput x; @rput y; @rput c; @rput w; @rput B; @rput sigma_sq
	@suppress R"""
		out = optrdd(X=x, W=w, Y=y, estimation.point=c,
					max.second.derivative=B, sigma.sq=sigma_sq,
					optimizer="quadprog", verbose=FALSE)
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
function compute_opt_bw(data::RDData)
	x = data.x; y = data.y; c = 0
	@rput x; @rput y; @rput c
	R"""
		bw_selection = rdrobust::rdbwselect_2014(y, x, c=c, bwselect="IK")
		opt_bw = bw_selection[["bws"]]
	"""
	@rget opt_bw
	return opt_bw[1]
end

function llrrdd(data::RDData, kernel, bw, ret_weights=false)
	x = data.x; y = data.y; c=0
	n = length(data.x)
	@rput y; @rput x; @rput c; @rput bw
	R"""
		model <- rdd::RDestimate(y~x, data.frame(y=y, x=x), cutpoint=c, bw=bw)
		est <- model$est[1]
	"""
	@rget est
	if ret_weights
		Δ = data.x .- c; w = data.x .>= c
		X = [ones(n) w min.(Δ, 0) max.(Δ, 0)]
		K = Diagonal(kernel.(abs.(Δ)./bw))
		H = inv(X'*K*X)*X'*K
		weights = H[2, :]
		return est, weights
	end
	return est
end
