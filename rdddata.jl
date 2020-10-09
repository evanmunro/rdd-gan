using DataFrames, GLM, Statistics

struct RDData
     y::Array{Float64, 1}
     x::Array{Float64, 1}
     sig2::Float64
     d0::Int
     d1::Int
     xx0::Array{Float64, 1}
     xx1::Array{Float64, 1}
     mu0::Array{Float64, 1}
     mu1::Array{Float64, 1}
     n0::Array{Int, 1}
     n1::Array{Int, 1}
end

function RDData(y, x)
      xx::Array{Float64, 1} = sort!(unique(x))
      xx0 = xx[xx .< 0]
      xx1 = xx[xx .> 0]
      d0 = length(xx0)
      d1 = length(xx1)
      data = DataFrame(y=y, x=x, w= x.>0)
      yhat = predict(lm(@formula(y~x*w),data))
      sig2 = mean((data.y - yhat).^2) * length(data.x) / (length(data.x) - 4)
      mu = zeros(Float64, length(xx))
      n = zeros(Int, length(xx))
      for i in 1:length(xx)
          mu[i] = mean(y[x.==xx[i]])
          n[i] = length(x[x.==xx[i]])
      end

      mu0 = mu[xx.<0]
      mu1 = mu[xx.>0]
      n0 = n[xx.<0]
      n1 = n[xx.>0]
      return RDData(y, x, sig2, d0, d1, xx0, xx1, mu0, mu1, n0, n1)
 end

function get_gamma(rd::RDData, weights)
      γ0 = zeros(rd.d0)
      γ1 = zeros(rd.d1)
      for i in 1:rd.d0
            γ0[i] = first(weights[data.x.== data.xx0[i]])
      end
      for i in 1:rd.d1
            γ1[i] = first(weights[data.x .== data.xx1[i]])
      end
      return γ0, γ1
end
