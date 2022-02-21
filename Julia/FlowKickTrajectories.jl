function flow_trajectory(model::Function, p, τ, starting_values)
    tspan = (0.0, τ)
    tvals = 0.0:0.1:τ
    prob = ODEProblem(model, starting_values, tspan, p)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    return solend = sol(tvals)
end

function flow_kick_trajectory(model::Function, p, k, τ, N, starting_values)
    tspan = (0.0, τ*N)
    tvals = 0.0:τ/50:τ*N

    function kicks(integrator)
        integrator.u[1] = integrator.u[1]+k
    end

    cb = PeriodicCallback(kicks, τ, initial_affect = false)
    prob = ODEProblem(model, starting_values, tspan, p)
    sol = DifferentialEquations.solve(prob, callback = cb,  reltol = 1e-8)
    return solend = sol(tvals)
end