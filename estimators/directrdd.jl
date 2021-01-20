using JuMP, Ipopt, Statistics, Optim, GLM, Distributions, Random

struct cvData
   ydot
   xdot
   ndot
   ys
   xs
   ns
end

function cvData(y, x; p=0.15)
    x = x ./ maximum(x)

    x̄ = unique(x)
    xcutoff = quantile(x̄, p)
    xdot = x̄[x̄.< xcutoff]
    n̄ = zeros(length(x̄))
    ȳ = zeros(length(x̄))
    for i in 1:length(x̄)
        ȳ[i] = mean(y[ x.== x̄[i]])
        n̄[i] = sum(x .== x̄[i])
    end

    xdot = x̄[x̄.<=xcutoff]
    ndot = n̄[x̄.<=xcutoff]
    ydot = ȳ[x̄.<=xcutoff]
    xs = [ x̄[x̄ .> xdot[1]] ]
    ns = [ n̄[x̄.==xdot[1]] ]
    ys = [ ȳ[x̄ .> xdot[1]] ]

    for i in 2:length(xdot)
        push!(xs,  x̄[x̄ .> xdot[i]] )
        push!(ns, n̄[x̄ .> xdot[i]] )
        push!(ys, ȳ[x̄ .> xdot[i]] )
    end
    return cvData(ydot, xdot, ndot, ys, xs, ns)
end

function rddCV(y, x, γ)
    y1 = y[x .> 0]; x1 = x[x .>0]
    y0 = y[x .< 0]; x0 = x[x.< 0]
    mu10 = predictCutoffCV(y1, x1, γ)
    mu00 = predictCutoffCV(y0, -1 .*x0, γ)
    return mu10 - mu00
end
function rddCV(y, x)
    y1 = y[x .> 0]; x1 = x[x .>0]
    y0 = y[x .< 0]; x0 = x[x.< 0]
    mu10 = predictCutoffCV(y1, x1)
    mu00 = predictCutoffCV(y0, -1 .*x0)
    return mu10 - mu00
end

function cvBWSelect(y, x)
    f = γ -> rddCVObjective(γ, y, x)
    res = Optim.optimize(f, 0.05, 0.9; rel_tol=1e-3)
    return res.minimizer
end


function rddCVSingle(y, x)
    return rddCV(y, x, cvBWSelect(y, x))
end

function rbwRDD(y, x)
    y1 = y[x .> 0]; x1 = x[x .>0]
    y0 = y[x .< 0]; x0 = x[x.< 0]
    mu10 = predictRand(y1, x1)
    mu00 = predictRand(y0, -1 .*x0)
    return mu10 - mu00
end

function predictRand(y, x)
    γ = rand(Uniform(0.05, 0.6), 1)[1]
    return triangularWLS(y, x, ones(length(x)), γ, 0.0)
end

function predictCutoffCV(y, x, γ)
    return triangularWLS(y, x, ones(length(x)), γ, 0.0)
end


function predictCutoffCV(y, x)
    data = cvData(y, x)
    f = γ -> rddCVObjective(γ, data)
    res = Optim.optimize(f, 0.05, 0.9; rel_tol=1e-3)
    γ = res.minimizer
    println(γ)
    return triangularWLS(y, x, ones(length(x)), γ, 0.0)
end

function triangularWLS(y, x, n, γ, x_eval)
    Δx = (x .- x_eval)
    ω = (1 .- Δx./ γ) .* (Δx .< γ).*n
    model = glm(@formula(y ~ x ), DataFrame(y=y, x=x), Normal(), IdentityLink(), wts=ω)
    return GLM.predict(model, DataFrame(x=[x_eval]))[1]
end

function rddCVObjective(γ, y, x)
    y1 = y[x .> 0]; x1 = x[x .>0]
    y0 = y[x .< 0]; x0 = x[x.< 0]
    data1 = cvData(y1, x1)
    data0 = cvData(y0, -1 .*x0)
    return (length(y1)*rddCVObjective(γ, data1) + length(y0)*rddCVObjective(γ, data0))/(length(y))
end

function rddCVObjective(γ, data::cvData)
    obj = zeros(length(data.xdot))
    for i in 1:length(data.xdot)
        yhat = triangularWLS(data.ys[i], data.xs[i], data.ns[i], γ, data.xdot[i])
        obj[i] = data.ndot[i]*(data.ydot[i] - yhat)^2
    end
    return sum(obj)/sum(data.ndot)
end
