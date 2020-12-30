using JuMP, Ipopt, Statistics, GLMNet, Optim

function directRDD(y, x, w)
    return 0
end

function localData(y, x)
    xcutoff = quantile(x, 0.05)
    ystar = y[ x .<= xcutoff]
    xstar = x[ x .<= xcutoff]
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

function dWeightsPoly(y, x; degree=1)
    f = θ -> objectivepoly(θ, y, x)
    res =Optim.optimize(f, [0.001, -0.1])
    println("obj min: ", res.minimum)
    θ = res.minimizer
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
    ystar, ys, xs = localData(y, x)
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

    return value(β₀), value(h), sum(value(β₀)*y[ x.<= value(h)])
end
