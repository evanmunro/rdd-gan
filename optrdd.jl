
using JuMP, MosekTools

include("estimators.jl")
include("rdddata.jl")

function optrdd(x, y, num_buckets=10000.0)
    B = estimate_B2(x, y)
    xdiscr = DiscretizedRunningVariable(x, num_buckets)
    X = support(xdiscr); h = xdiscr.h; d = length(X); n= xdiscr.weights; ixc = xdiscr.ixc
    c = 0.0; W = X .> c
    σ2 = zeros(d) .+ estimate_σ2(x ,y)
    model =  Model(Mosek.Optimizer)

    @variable(model, G[1:d])
    @variable(model, f[1:d])
    #@variable(model, f0[1:d])
    #@variable(model, f1[1:d])
    @variable(model, l1 >=0)
    @variable(model, l2)
    @variable(model, l3)
    @variable(model, l4)
    @variable(model, l5)

    for i in 1:d
        @constraint(model, G[i] == 2*B*f[i] + l2*(1-W[i]) + l3*W[i] + l4*(X[i] - c) + l5*(W[i] - 0.5)*(X[i]-c))
    end

    #for i in 1:d
        #@constraint(model, f[i] == (1-W[i])*f0[i] + W[i]*f1[i])
    #end

    @constraint(model, f[ixc] ==0)
    @constraint(model, (f[ixc] - f[ixc-1])/h[ixc] == 0)
    @constraint(model, (f[ixc+1] - f[ixc])/h[ixc] == 0)

    #for i in 2:d
    #    @constraint(model, (f1[i] - f1[i-1])/h - l1*2 <=0)
    #    @constraint(model, (f1[i] - f1[i-1])/h  >= 0)
    #    @constraint(model, (f0[i] - f0[i-1])/h - l1*2 <=0)
    #    @constraint(model, (f0[i] - f0[i-1])/h  >= 0)
    #end

    for i in 3:d
        @constraint(model, (f[i] - 2*f[i-1] + f[i-2])/(h[i-1]*h[i-2]) - l1 <=0)
        @constraint(model, (f[i] - 2*f[i-1] + f[i-2])/(h[i-1]*h[i-2]) + l1 >=0)
        #@constraint(model, (f0[i] - 2*f0[i-1] + f0[i-2])/h^2 - l1 <=0)
        #@constraint(model, (f0[i] - 2*f0[i-1] + f0[i-2])/h^2 + l1 >=0)
    end

    @objective(model, Min, 1/4*(sum(n[i]*G[i]^2/σ2[i] for i in 1:d)) + l1^2 - l2 + l3)

    optimize!(model)
    println("l1: ", value(l1))
    println("l2: ",  value(l2))
    println("l3: ", value(l3))
    println("l4: ", value(l4))
    println("l5: ", value(l5))
    println("obj: ", objective_value(model))

    γ_xx = -value.(G)./(2 .*σ2)
    γ = γ_xx[xdiscr.xmap]
    return (τ = sum(γ.*y), γ = γ)
end
