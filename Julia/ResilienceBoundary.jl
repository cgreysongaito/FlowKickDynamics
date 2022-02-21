#For each harvest (kick) size, the script approximates the minimum time it takes to population to recover by that magnitude.
function res_bound(model::Function, p, krange, stablept, unstablept, x0stepsize)
    #Create a range of kick values 
    #Create an array to store points on the resilience boundary.
    ResilienceCurveτ = zeros(length(krange))
    for (ki, knum) in enumerate(krange)
        if unstablept < stablept
            x0range = unstablept+x0stepsize-knum:x0stepsize:stablept-x0stepsize #identify a set of x values at which to start flowing (min computed over this set)
        else
            x0range = stablept+x0stepsize:x0stepsize:unstablept-x0stepsize-knum
        end
        TimeToFlow = zeros(length(x0range))
            for (x0i, x0num) in enumerate(x0range)
                TimeToFlow[x0i] = quadgk(x -> 1/model(x, p), x0num+knum, x0num)[1] #calculate time to flow from x0+k to x0 #USES the inverse of the model
            end
        ResilienceCurveτ[ki] = minimum(TimeToFlow)
    end
    return [ResilienceCurveτ, abs.(collect(krange))]
end