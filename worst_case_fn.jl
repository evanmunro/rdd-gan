using JuMP, Ipopt
include("rdddata.jl")
function worst_case_mu(rd::RDData, γ0, γ1, M)
    model =  Model(Gurobi.Optimizer)
    @variable(model, mu0[1:rd.d0])
    @variable(model, mu1[1:rd.d1])
    mu00=0
    mu10=0
    #@variable(model, mu0[1:rd.d0])
    #@variable(model,mu1[1:rd.d1])

    #@variable(model, mu00 )
    #@variable(model,  mu10)


    # second derivative constraints at cutoff
    @constraint(model, -M <= (mu1[2] - 2*mu1[1] + mu10 )/(rd.xx1[1]*(rd.xx1[2]
                                                            -rd.xx1[1])) <= M )
    @constraint(model, -M <= (mu00 - 2*mu0[rd.d0] + mu0[rd.d0-1])/(rd.xx0[rd.d0]*(rd.xx0[rd.d0]
                                                            -rd.xx1[rd.d0-1])) <= M )

    #second derivative constraints away from the cutoff
    for i in 3:rd.d1
        @constraint(model, -M <= (mu1[i] - 2*mu1[i-1] + mu1[i-2] )/((rd.xx1[i] - rd.xx1[i-1])*(rd.xx1[i-1]
                                                                -rd.xx1[i-2])) <= M )
    end

    for i in 1:(rd.d0-2)
        @constraint(model, -M <= (mu0[i+2] - 2*mu0[i+1] + mu0[i])/((rd.xx0[i+2] - rd.xx0[i+1])*(rd.xx0[i+1]
                                                                -rd.xx0[i])) <= M )
    end

    #first derivative constraints at the cutoff
    @constraint(model, (mu1[1] - mu10)/rd.xx1[1] == 0)
    @constraint(model, (mu00 - mu0[rd.d0])/rd.xx0[rd.d0] == 0)

    # first derivative constraints away from the cutoff

    #Bounding box constraints away from the cutoff

    # maximize the absolute bias
    @objective(model, Max, (sum(rd.n0[i]*γ0[i]*mu0[i] for i in 1:rd.d0) +
                            sum(rd.n1[i]*γ1[i]*mu1[i] for i in 1:rd.d1)
                            - mu10 + mu00))

    optimize!(model)
    max_bias = objective_value(model)
    max_mu = [value.(mu0)..., value.(mu1)...]

    @objective(model, Min, (sum(rd.n0[i]*γ0[i]*mu0[i] for i in 1:rd.d0) +
                            sum(rd.n1[i]*γ1[i]*mu1[i] for i in 1:rd.d1)
                            - mu10 + mu00))

    optimize!(model)
    min_bias = objective_value(model)
    min_mu = [value.(mu0)..., value.(mu1)...]

    max_abs_bias = max(abs(max_bias), abs(min_bias))

    return (max_abs_bias = max_abs_bias,
            max_bias = max_bias,
            min_bias = min_bias,
            max_mu = max_mu,
            min_mu = min_mu)
end
