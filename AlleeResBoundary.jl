#THIS SCRIPT PLOTS THE FLOW-KICK RESILIENCE BOUNDARY FOR ALLEE POPULATIONS 1 AND 2
#For each harvest (kick) size, the script approximates the minimum time it takes to population to recover by that magnitude.

include("JuliaPackages.jl")

@with_kw mutable struct AlleePar
    K=100.0  #carrying capacity
    A=20.0 #critical Allee threshold
    α=0.0 #first parameter in quadratic modifier
    β=0.0 #second parameter in quadratic modifier
    γ=1.0 #third parameter in quadratic modifier - set at 1.0 for when don't want quadratic modifier
end

function Allee_model(x, p)
    @unpack K, A, α, β, γ = p
    return (x * (1 - (x/K)) * ((x/A) - 1) * ((α * x^2) - (β * x) + γ))
end

function res_bound(model::Function, p, krange, stable, unstable, x0stepsize)
    #Create a range of kick values 
    #Create an array to store points on the resilience boundary.
    ResilienceCurveτ = zeros(length(krange))
    for (ki, knum) in enumerate(krange)
        if unstable < stable
            x0range = unstable+x0stepsize-knum:x0stepsize:stable-x0stepsize #identify a set of x values at which to start flowing (min computed over this set)
        else
            x0range = stable+x0stepsize:x0stepsize:unstable-x0stepsize-knum
        end
        TimeToFlow = zeros(length(x0range))
            for (x0i, x0num) in enumerate(x0range)
                TimeToFlow[x0i] = quadgk(x -> 1/model(x, p), x0num+knum, x0num)[1] #calculate time to flow from x0+k to x0 #USES the inverse of the model
            end
        ResilienceCurveτ[ki] = minimum(TimeToFlow)
    end
    return [ResilienceCurveτ, abs.(collect(krange))]
end

let
    data1 = res_bound(Allee_model, AlleePar(), -79.8:0.1:0.0, 100.0, 20.0, 0.1)
    data2 = res_bound(Allee_model,AlleePar(α = 0.0002, β=0.024, γ=1.4), -79.8:0.1:0.0, 100.0, 20.0, 0.1)
    test = figure()
    plot(data1[1],data1[2])
    plot(data2[1],data2[2])
    xlabel("flow time (τ)")
    ylabel("kick size (κ)")
    xlim(0,5)
    ylim(0,80)
    return test
end
