using JuMP, Ipopt, Statistics

function directRDD(y, x, w)
    return 0
end

function localData(y, x)
    xcutoff = quantile!(x, 0.25)
    ystar = y[ x .<= xcutoff]
    xstar = x[ x .<= xcutoff]
    xs = [ x[ x .< xlim] for xlim in xstar ]
    ys = [ y[ x .< xlim] for xlim in xstar ]
    return ystar, ys, xs
end


function directWeightsConstant(y, x)
    ystar, ys, xs = localData(y, x)
    d = length(ystar)
    model = Model(Ipopt.Optimizer)
    @variable(model, β₀ >= 0)
    @variable(model, h >= 0, start=1)
    #@NLconstraint(model, β₀*sum(x[i] <= h for i in length(x)) == 1)

    @variable(model, obj[1:d])

    for i in 1:d
        @NLconstraint(model, obj[i] == ystar[i] -
            sum(β₀*ys[i][j]*(xs[i][j]<=h) for j in 1:length(ys[i]))/sum(xs[i][j]<=h for j in 1:length(ys[i])))
    end

    @objective(model, Min, sum(obj[i]^2 for i in 1:d))
    optimize!(model)

    return value(β₀), value(h), sum(value(β₀)*y[ x.<= value(h)])
end

function directWeightsLinear(y, x)
    return 0
end

function directWeightsPoly(y, x)
    return 0
end
