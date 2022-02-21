#THIS SCRIPT APPROXIMATES THE AREA ABOVE A FLOW-KICK RESILIENCE BOUNDARY AND BELOW ITS HORIZONTAL ASYMPTOTE.
#For each k, find the minimum time it takes to flow forward that distance

function res_area_prep(model::Function, p, krange, stablept)
    ResilienceCurve_τ=zeros(length(krange))
    if minimum(krange) < 0.0
        xmin=[stablept-step(krange)]
    else
        xmin=[stablept+step(krange)]
    end
    for (ki, knum) in enumerate(krange)
       append!(xmin, find_zero(x->model(x+knum,p)-model(x,p),xmin[ki]))
   #find root eqn. fzero function needs the function and initial guess
#     initial guess should be between 20 and 100. start with x0=99.9 and
#     then use previous answer for next x0
        ResilienceCurve_τ[ki] = quadgk(x -> 1/model(x,p),xmin[ki+1]+knum,xmin[ki+1])[1] #the tau that goes with k
    end
    return [ResilienceCurve_τ, abs.(collect(krange))]
end

function res_area(res_curve_data, basin)
    τ = reverse(res_curve_data[1])
    k = reverse(res_curve_data[2])
    area1=trapz(τ,basin .- k) #trap rule d(tau) with unequal intervals
    area2=trapz(k,τ) #trap rule d(kick) with equal intervals
    return (area1+area2)/2
end