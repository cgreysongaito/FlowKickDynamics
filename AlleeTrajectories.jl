#THIS SCRIPT CALCULATES FLOW-KICK TRAJECTORIES FOR A HARVESTED POPULATION
include("JuliaPackages.jl")

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

function allee!(du, u, p, t)
    du[1] = u[1] * (1.0 - (u[1]/100.0)) * ((u[1]/20.0) - 1.0)
end

let 
    data1 = flow_kick_trajectory(allee!,-12.0, 0.25, 100.0, 64)
    data2 = flow_kick_trajectory(allee!,-40.0, 10/12, 100.0, 15)
    flowkicktrajectory = figure()
    plot(data1.t, data1.u)
    plot(data2.t, data2.u)
    xlabel("Time (years)")
    ylabel("Population Size \n (kt biomass)")
    return flowkicktrajectory
end