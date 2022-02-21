function newton(coupledvar::Function, p, k, τ, x0)
    #Set up a place to hold Newton iterates x_n 
    Newton_pts=[x0]
    converge = 1

    for i in 1:10
        x = Newton_pts[i]
        prob = ODEProblem(coupledvar, [x, 1.0], (0.0,τ), p)
        sol = DifferentialEquations.solve(prob,  reltol = 1e-8)
        post_flow = sol[end][1] #this is phi^tau of x
        Vsoln = sol[end][2] #this is the derivative of phi^tau with respect to x
        
        DF = Vsoln-1 #Compute DF
        F_of_x = post_flow + k - x #Compute F(x)
        newx = x - F_of_x/DF
        append!(Newton_pts, [newx])
        if abs(F_of_x) < 10^(-10) #check whether x is nearly at equilibrium
            converge = 0
            break
        end
    end
    return [Newton_pts[end], converge]
end

#This script traces out a branch of fixed points for a flow-kick map using Newton's method and parameter continuation

#use equilibria from deterministic model to set as guesses then run for loop over each guess to retreive each branch

# NOTE krange must always be in sequential order (negative numbers before 0.0 then positive numbers)
function branch(coupledvar::Function, p, krange, τ, initguess)
    fxpdts = [initguess]
    if minimum(krange) < 0.0
        orderedkrange = reverse(krange)
    else
        orderedkrange = krange
    end
    finalki=0
    for (ki, knum) in enumerate(orderedkrange)
        fixed_point = newton(coupledvar, p, knum, τ, fxpdts[end])
        if fixed_point[2] == 1.0
            finalki = ki-1
            break
        else
            finalki = ki
            append!(fxpdts, fixed_point[1])
        end
    end
    return [abs.(orderedkrange[1:finalki]),fxpdts[2:end]]
end