function flow_trajectory(model::Function, p, τ, x0, y0)
    tspan = (0.0, τ)
    tvals = 0.0:0.1:τ
    prob = ODEProblem(model, [x0, y0], tspan, p)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    return solend = sol(tvals)
end

function flow_kick_trajectory(model::Function, p, k, τ, x0, N)
    tspan = (0.0, τ*N)
    tvals = 0.0:τ/4:τ*N

    function kicks(integrator)
        integrator.u[1] = integrator.u[1]+k
    end

    cb = PeriodicCallback(kicks, τ, initial_affect = false)
    prob = ODEProblem(model, x0, tspan, p)
    sol = DifferentialEquations.solve(prob, callback = cb,  reltol = 1e-8)
    return solend = sol(tvals)
end