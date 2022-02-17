#This script traces out a branch of fixed points for a flow-kick map using Newton's method and parameter continuation

#use equilibria from deterministic model to set as guesses then run for loop over each guess to retreive each branch

function guess_fp(fxdpts, klist, kstep)
    if fxdpts[end] == fxdpts[end-1]
        return fxdpts[end]
    else
        return fxdpts[end] + kstep * (fxdpts[end] - fxdpts[end-1])/(klist[end] - klist[end-1])
    end
end

function branch(kmin, kmax, kstep, τ, initguess)
    krange = -kmin:-kstep:-kmax #what about negative or positive kicks?
    fxpdts = [initguess, initguess]
    klist = [kmin, kmin]
    for knum in krange
        fixed_point = newton(knum, τ, guess_fp(fxpdts, klist, kstep))
        if fixed_point[2] == 1.0
            break
        else
            append!(fxpdts, fixed_point[1])
            append!(klist, -knum)
        end
    end
    return [klist[3:end],fxpdts[3:end]]
end

