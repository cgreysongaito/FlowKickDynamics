function CoupledVar!(du, u, p, t)
    #gives ODE for fish population 1, along with varational equation to be solved simultaneously
    du[1] = u[1] * (1.0 - (u[1]/100.0)) * ((u[1]/20.0) - 1.0)
    du[2] = (-1+3*(u[1]/25) - 3*(u[1]^2/2000))*u[2]
end

function newton(k, τ, x0)
    #Set up a place to hold Newton iterates x_n 
    Newton_pts=[x0]
    converge = 1

    for i in 1:10
        x = Newton_pts[i]
        prob = ODEProblem(CoupledVar!, [x, 1.0], (0.0,τ))
        sol = DifferentialEquations.solve(prob,  reltol = 1e-8)
        post_flow = sol[end][1] #this is phi^tau of x
        Vsoln = sol[end][2] #this is the derivative of phi^tau with respect to x
        
        DF = Vsoln-1 #Compute DF
        F_of_x = post_flow + k - x #Compute F(x)
        newx = x - F_of_x/DF #x=x-inv(DF)*F_of_x
        append!(Newton_pts, [newx])
    
        if abs(F_of_x) < 10^(-10) #check whether x is nearly at equilibrium
            converge = 0
            break
        end
    end
    return [Newton_pts[end], converge]
end