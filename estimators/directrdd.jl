using JuMP, Ipopt, Statistics, Optim, GLM, Distributions, Random

struct cvData
   ydot
   xdot
   ndot
   ys
   xs
   ns
end

function cvData(y, x; p=0.05)
    x = x ./ maximum(x)

    x̄ = unique(x)
    xcutoff = quantile(x̄, p)
    xdot = x̄[x̄.< xcutoff]
    counts = countmap(x)
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
function rddCV(y, x)
    y1 = y[x .> 0]; x1 = x[x .>0]
    y0 = y[x .< 0]; x0 = x[x.< 0]
    mu10 = predictCutoffCV(y1, x1)
    mu00 = predictCutoffCV(y0, -1 .*x0)
    return mu10 - mu00
end

function rbwRDD(y, x)
    y1 = y[x .> 0]; x1 = x[x .>0]
    y0 = y[x .< 0]; x0 = x[x.< 0]
    mu10 = predictRand(y1, x1)
    mu00 = predictRand(y0, -1 .*x0)
    return mu10 - mu00
end

function predictRand(y, x)
    γ = rand(Uniform(0.17, 0.99), 1)[1]
    return triangularWLS(y, x, ones(length(x)), γ, 0.0)
end


function predictCutoffCV(y, x)
    data = cvData(y, x)
    f = γ -> rddCVObjective(γ, data)
    res = Optim.optimize(f, 0.1, 1.0; rel_tol=1e-3)
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

function rddCVObjective(γ, data::cvData)
    obj = zeros(length(data.xdot))
    for i in 1:length(data.xdot)
        yhat = triangularWLS(data.ys[i], data.xs[i], data.ns[i], γ, data.xdot[i])
        obj[i] = data.ndot[i]*(data.ydot[i] - yhat).^2
    end
    return sum(obj)
end

function localData(y, x)
    ystar = y[ x .<= xcutoff]
    xs = [ x[ x .> xlim] for xlim in xstar ]
    ys = [ y[ x .> xlim] for xlim in xstar ]
    return ystar, xstar, ys, xs
end

function rddLinear(y, x)
    return 0
end

function rddLasso(y, x)
    y1 = y[x .> 0]; x1 = x[x .>0]
    y0 = y[x .< 0]; x0 = x[x.< 0]
    h1 = selectWindow(y1, x1)
    h0 = selectWindow(y0, x0)
    mu10 = lassoFit(y1[abs.(x1) .< h1], x1[abs.(x1) .< h1], 0.0)
    mu00 = lassoFit(y0[abs.(x0) .< h0], x0[abs.(x0) .< h0], 0.0)
    return mu10 - mu00

end

function windowObjective(h, y, x)
    ystar, xstar, ys, xs = localData(y, x)
    d = length(ystar)
    obj = zeros(d)
    for i in 1:d
        inWindow = abs.(xs[i] .- xstar[i]) .<= h
        obj[i] = ystar[i] - lassoFit(ys[i][inWindow], xs[i][inWindow], xstar[i])
    end

    return sum(obj.^2)
end

function selectWindow(y, x)
    f = h -> windowObjective(h, y, x)
    res = Optim.optimize(f, 0.05, 0.5; rel_tol=1e-3)
    h = res.minimizer
end

function lassoFit(y, x, xeval)
    X = [x x.^2 x.^3]
    cv = glmnetcv(X, y; alpha=0.1)
    return GLMNet.predict(cv, [xeval xeval.^2 xeval.^3])[1]
end

function rddConstant(y, x)
    mu1 = dWeightsConstant(y[ x.>0], x[x.>0])
    mu0 = dWeightsConstant(y[x.<0], -1 .*x[x.<0])
    return mu1 - mu0
end

function rddPoly(y, x)
    mu1 = dWeightsPoly(y[ x.>0], x[x.>0])
    mu0 = dWeightsPoly(y[x.<0], -1 .*x[x.<0])
    println(mu1)
    println(mu0)
    return mu1 - mu0
end

function rddPoly2(y, x)
    a, mu1 = optimalWeights(y[ x.>0], x[x.>0], 0.005, 2)
    a, mu0 = optimalWeights(y[x.<0], -1 .*x[x.<0], 0.005, 2)
    println(mu1)
    println(mu0)
    return mu1 - mu0
end

function optimalWeights(y, x, h, degree)
    ystar, xstar, ys, xs = localData(y, x)
    model = Model(Ipopt.Optimizer)
    @variable(model, α)
    @variable(model, β₀, start=1)
    @variable(model, β[1:degree])
    @constraint(model, α==0)
    #@NLconstraint(model, β₀*sum(x[i] <= h for i in 1:length(x)) == 1)

    @variable(model, obj[1:length(ystar)])

    for i in 1:length(ystar)
        @NLconstraint(model, obj[i] == ystar[i] - α -
            sum((β₀ + sum(β[p]*(xs[i][j] - xstar[i])^p for p in 1:degree))*((xs[i][j] - xstar[i]) >h)*ys[i][j]
            for j in 1:length(ys[i]))/sum((β₀ + sum(β[p]*(xs[i][j] - xstar[i])^p for p in 1:degree))*((xs[i][j] - xstar[i]) >h)
                                                for j in 1:length(ys[i])))
    end

    @objective(model, Min, sum(obj[i]^2 for i in 1:length(ystar))/length(ystar))
    optimize!(model)

    βv = value.(β)
    wx = (value(β₀) .+ sum(βv[p].* x.^p for p in 1:degree)).* (x.< h)
    μ̂ = sum(wx.*y)/sum(wx) + value(α)
    println(value(β₀))
    println(βv)
    return  objective_value(model), μ̂
end

function dWeightsPoly(y, x; degree=1)
    f = h -> objectivepoly(h, y, x)
    res =Optim.optimize(f, 0.05, 0.5; rel_tol=1e-3)
    println("obj min: ", res.minimum)
    h = res.minimizer
    println(θ)
    h = θ[1]
    w(x) = sum(θ[i]*x^(i-1) for i in 2:length(θ))
    β0 = 1/sum(x .<= h) - sum(w.(x[x.<=h]))/sum(x .<= h)
    println(sum(β0 .+ w.(x[x.<=h])))
    return sum((β0 .+ w.(x[x.<=h])).*y[ x.<= h])
end

function objectivepoly(θ, y, x)
    ystar, xstar, ys, xs = localData(y, x)
    d = length(ystar)
    h = θ[1]
    w(x) = sum(θ[i]*x^(i-1) for i in 2:length(θ))
    β0 = 1/sum(x .<= h) - sum(w.(x[x.<=h]))/sum(x .<= h)
    obj = zeros(d)
    for i in 1:d
        obj[i] = ystar[i] -
        sum(ys[i][j]*(β0 +
           w(xs[i][j]))*(xs[i][j]<=h) for j in 1:length(ys[i]))/sum((β0 .+ w.(xs[i])).*(xs[i].<=h))
    end

    return sum(obj.^2)
end

function dWeightsConstant(y, x)
    f = h -> objective(h, y, x)
    res = Optim.optimize(f, minimum(x), maximum(x); rel_tol=1e-3)
    h = res.minimizer
    β = 1/sum(x .<= h)
    return sum(β*y[ x.<= h])
end

function objective(h, y, x)
    ystar, xstar, ys, xs = localData(y, x)
    β = 1/sum(x .<= h)
    d = length(ystar)
    obj = zeros(d)
    for i in 1:d
        obj[i] = ystar[i] -
        sum(β*ys[i][j]*(xs[i][j]<=h) for j in 1:length(ys[i]))/sum(xs[i].<=h)
    end

    return sum(obj.^2)
end


function directWeightsPoly(y, x)

    d = length(ystar)
    p = 3
    model = Model(Ipopt.Optimizer)
    @variable(model, β₀)
    @variable(model, β[1:p])
    #@NLconstraint(model, β₀*sum(x[i] <= h for i in 1:length(x)) == 1)

    @variable(model, obj[1:d])

    for i in 1:d
        @NLconstraint(model, obj[i] == ystar[i] -
            sum(β₀*ys[i][j]*(xs[i][j]<=h) for j in 1:length(ys[i]))/sum(xs[i][j]<=h for j in 1:length(ys[i])))
    end

    @objective(model, Min, sum(obj[i]^2 for i in 1:d))
    optimize!(model)

    return objective_value(model)
end
