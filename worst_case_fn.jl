using JuMP, Ipopt
import RegressionDiscontinuity


#unction add_value_bounds(model, variable, y)
#end

#function add_monotonicity_bounds(model, variable, M1, y)
#end

#function add_second_derivative_bounds(model, variable, M1)

#end

function worst_case_mu(rd::RegressionDiscontinuity.RDData, γ0, γ1, M1, M2, first_deriv=true)
    model =  Model(Gurobi.Optimizer)
    @variable(model, mu0[1:rd.d0])
    @variable(model, mu1[1:rd.d1])
    mu00=0.0
    mu10=0.0
    #@variable(model, mu0[1:rd.d0])
    #@variable(model,mu1[1:rd.d1])

    #@variable(model, mu00 )
    #@variable(model,  mu10)


    # second derivative constraints at cutoff
    @constraint(model, -M2 <= (mu1[2] - 2*mu1[1] + mu10 )/(rd.xx1[1]*(rd.xx1[2]
                                                            -rd.xx1[1])) <= M2 )
    @constraint(model, -M2 <= (mu00 - 2*mu0[rd.d0] + mu0[rd.d0-1])/(rd.xx0[rd.d0]*(rd.xx0[rd.d0]
                                                            -rd.xx1[rd.d0-1])) <= M2 )

    #second derivative constraint across cutoff
    @constraint(model, -M2 <= (mu1[1] - 2*mu10 + mu0[rd.d0])/(rd.xx1[1]*rd.xx0[rd.d0]) <= M2 )


    #second derivative constraints away from the cutoff
    for i in 3:rd.d1
        @constraint(model, -M2 <= (mu1[i] - 2*mu1[i-1] + mu1[i-2] )/((rd.xx1[i] - rd.xx1[i-1])*(rd.xx1[i-1]
                                                                -rd.xx1[i-2])) <= M2 )
    end

    for i in 1:(rd.d0-2)
        @constraint(model, -M2 <= (mu0[i+2] - 2*mu0[i+1] + mu0[i])/((rd.xx0[i+2] - rd.xx0[i+1])*(rd.xx0[i+1]
                                                                -rd.xx0[i])) <= M2 )
    end

    #first derivative constraints at the cutoff
    @constraint(model, (mu1[1] - mu10)/rd.xx1[1] == 0)
    @constraint(model, (mu00 - mu0[rd.d0])/rd.xx0[rd.d0] == 0)

    # first derivative constraints away from the cutoff

    if first_deriv
        for i in 2:rd.d1
            @constraint(model, 0 <= (mu1[i]- mu1[i-1])/(rd.xx1[i]-rd.xx1[i-1]) <= M1)
        end

        for i in 1:(rd.d0-1)
            @constraint(model, 0 <= (mu0[i+1] - mu0[i])/(rd.xx0[i+1]-rd.xx0[i]) <= M1 )
        end
    end
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
