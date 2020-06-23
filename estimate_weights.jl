function estimate_gamma(rd::RDData)
    convg = false
    alpha = 0
    mu0 = rd.mu[rd.c0-1]
    mu1 = rd.mu[rd.c1+1]
    println(mu0)
    println(mu1)
    gPrev = optimize_gamma(rd, mu0, mu1)
    while(!convg)
        println("step")
        mu0, mu1 = optimize_mu(rd, copy(gPrev))
        println(mu0)
        println(mu1)
        gNew =  optimize_gamma(rd, mu0, mu1)
        diff = mean(abs.(gPrev.-gNew))
        #println(diff)
        if diff < 0.00001
            convg = true
        end
        gPrev = alpha*gPrev + (1-alpha)*copy(gNew)
        #println(mean(abs.(gPrev.-gNew)))
    end
    return gPrev
end

function optimize_gamma(rd::RDData, mu0, mu1)
    model =  Model(with_optimizer(Ipopt.Optimizer))

    @variable(model, gamma[1:rd.d])

    @constraint(model,sum(rd.n[i]*gamma[i] for i in 1:rd.c0)==-1)
    @constraint(model,sum(rd.n[i]*gamma[i] for i in rd.c1:rd.d)==1)
    @constraint(model,sum(rd.n[i]*gamma[i]*rd.xx[i] for i in 1:rd.d)==0)

    @NLobjective(model,Min,
                    rd.sig2*sum(rd.n[i]*gamma[i]^2 for i in 1:rd.d)
                    + (sum(rd.n[i]*gamma[i]*rd.mu[i] for i in 1:rd.d)
                    - mu1+mu0)^2)
    optimize!(model)
    println(objective_value(model))
    return copy(value.(gamma))
end

function optimize_mu(rd::RDData, gamma)
    model =  Model(with_optimizer(Ipopt.Optimizer))

    @variable(model, mu0)
    @variable(model, mu1)

    # second derivative constraints
    @constraint(model,-rd.M <= (rd.mu[rd.c1+1] - 2*rd.mu[rd.c1] + mu1 )/rd.h^2 <= rd.M )
    @constraint(model,-rd.M <= (mu0 - 2*rd.mu[rd.c0] + rd.mu[rd.c0-1])/rd.h^2 <= rd.M )

    # first derivative constraints
    #@constraint(model, 0 <= (rd.mu[rd.c1] - mu1)/rd.h <= 1)
    #@constraint(model, 0 <= (mu0 - rd.mu[rd.c0])/rd.h <=1)

    @objective(model,Max,(sum(rd.n[i]*gamma[i]*rd.mu[i] for i in 1:rd.d)
                    - mu1+mu0)^2)

    optimize!(model)
    println(objective_value(model))
    return copy(value(mu0)), copy(value(mu1))

end




model = Model(with_optimizer(Ipopt.Optimizer))

@variable(model, gamma[1:d])
@variable(model, t)

#identification constraints (sum to 1)
@constraint(model, sum(n[i]*gamma[i] for i in 1:c_0)==1)
@constraint(model, sum(n[i]*gamma[i] for i in c_1:d)==-1)
@constraint(model, sum(n[i]*gamma[i]*values[i] for i in 1:d)==0)

for m in 1:M
    @constraint(t >= sig2*sum(n[i]*gamma[i]^2 for i in 1:d)
    + (sum(n[i]*gamma[i]*mu[i] for i in 1:d)-mu10[m]+mu00[m])^2) )
end
    

#second derivative constraint on either side of 0
#@constraint(model,-M <= (mu[c_1+1] - 2*mu[c_1] + mu10 )/h^2 <= M )
#@constraint(model, -M <= (mu00 - 2*mu[c_0] + mu[c_0-1])/h^2 <= M )

#first derivative constraint on either side of 0
#@constraint(model, 0 <= (mu[c_1] - mu10)/h <= 1)

#@constraint(model, 0 <= (mu00 - mu[c_0])/h <=1)

#constraint for max bias


#objective function
@NLobjective(model,Min,
    sig2*sum(n[i]*gamma[i]^2 for i in 1:d)
    + (sum(n[i]*gamma[i]*mu[i] for i in 1:d)-mu10+mu00)^2)

optimize!(model)
vg = value.(gamma)
