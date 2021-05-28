
function worst_case_mu(x, y, γ, γ0 = 0.0)
    model =  Model(Mosek.Optimizer)
    data  = DiscretizedRDData(y, x, 2000)
    bounds = AdaptiveFunctionClass(data)
    X = data.xx; d=length(X); h = data.h; n= data.weights; ixc = data.ixc
    @variable(model, mu0[1:d])
    @variable(model, mu1[1:d])
    c = 0.0; W = X .> c; j = data.xmap

    add_constraints!(model, mu0, bounds.μ₀, n.*(1 .- W), h)
    add_constraints!(model, mu1, bounds.μ₁, n.*W, h)
    ## add these constraints for model to have a solution
    #@constraint(model, mu0[ixc] == 0 )
    #@constraint(model, mu1[ixc] == 0 )
    #@constraint(model, (mu0[ixc+1] - mu0[ixc])/h[ixc] == 0)
    #@constraint(model, (mu1[ixc+1] - mu1[ixc])/h[ixc] == 0)

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

function add_second_derivative_constraint!(model, f, h, M2)
    d  = length(h) + 1
    for i in 3:d
        @constraint(model, M2 <= (f[i] - 2*f[i-1] + f[i-2])/(h[i-1]*h[i-2])<= M2)
    end
end


function bias_adjusted_gaussian_ci(se; maxbias=0.0, level=0.95)
    rel_bias = maxbias/se
    zz = fzero( z-> cdf(Normal(), rel_bias-z) + cdf(Normal(), -rel_bias-z) +  level -1,
        0, rel_bias - quantile(Normal(),(1- level)/2.1))
    zz*se
end
