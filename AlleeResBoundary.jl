#THIS SCRIPT PLOTS THE FLOW-KICK RESILIENCE BOUNDARY FOR ALLEE POPULATIONS 1 AND 2
#For each harvest (kick) size, the script approximates the minimum time it takes to population to recover by that magnitude.

include("JuliaPackages.jl")

#Create the inverse of the vector field function, i.e. dt/dx
function pop1(x)
   return 1/(x * (1 - (x/100)) * ((x/20) - 1))
end

function pop2(x)
    return 1/(x * (1 - (x/100)) * ((x/20) - 1) * ((0.0002 * x^2) - (0.024 * x) + 1.4))
end

function resboundAllee(popfun::Function)
    #Create a range of kick values 
    krange = -79.8:0.1:0;  #go close to distance to threshold, but not quite

    #Create an array to store points on the resilience boundary.
    ResilienceCurveτ = zeros(length(krange))
    ResilienceCurvek = zeros(length(krange))

    for (ki, knum) in enumerate(krange)
        x0range = 20.1-knum:0.1:99.9; #identify a set of x values at which to start flowing (min computed over this set)
        TimeToFlow = zeros(length(x0range))
            for (x0i, x0num) in enumerate(x0range)
                TimeToFlow[x0i] = quadgk(x -> popfun(x), x0num+knum, x0num)[1] #calculate time to flow from x0+k to x0
            end
        ResilienceCurveτ[ki] = minimum(TimeToFlow)
        ResilienceCurvek[ki] = abs(knum)
    end
    return [ResilienceCurveτ, ResilienceCurvek]
end

let
    data1 = resboundAllee(pop1)
    data2 = resboundAllee(pop2)
    test = figure()
    plot(data1[1],data1[2])
    plot(data2[1],data2[2])
    xlabel("flow time tau")
    ylabel("kick size |kappa|")
    xlim(0,5)
    ylim(0,80)
    return test
end