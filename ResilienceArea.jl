#THIS SCRIPT APPROXIMATES THE AREA ABOVE A FLOW-KICK RESILIENCE BOUNDARY AND BELOW ITS HORIZONTAL ASYMPTOTE.
include("JuliaPackages.jl")
#For each k, find the minimum time it takes to flow forward that distance
function pop1(x)
    return 1/(x * (1 - (x/100)) * ((x/20) - 1))
 end
function pop2(x)
    return 1/(x * (1 - (x/100)) * ((x/20) - 1) * ((0.0002 * x^2) - (0.024 * x) + 1.4))
end

function res_area_prep(model::Function)
    krange=-79.99:0.01:0; #go close to distance to threshold, but not quite
    ResilienceCurve_τ=zeros(length(krange))
    xmin=[99.99]
    for (ki, knum) in enumerate(krange)
    #der=@(x) -(x.*(1-(x/100)).*((x/20)-1))+((x+k).*(1-((x+k)/100)).*(((x+k)/20)-1));
    # der=@(x) -(x.*(1-(x/100)).*((x/20)-1).*(0.0002*x.^2-0.024*x+1.4))+((x+k).*(1-((x+k)/100)).*(((x+k)/20)-1).*(0.0002*(x+k).^2-0.024*(x+k)+1.4));
        append!(xmin, find_zero(x->1/model(x+knum)-1/model(x),xmin[ki]))
   #find root eqn. fzero function needs the function and initial guess
#     initial guess should be between 20 and 100. start with x0=99.9 and
#     then use previous answer for next x0
        ResilienceCurve_τ[ki] = quadgk(model,xmin[ki+1]+knum,xmin[ki+1])[1] #the tau that goes with k
    end
    return [ResilienceCurve_τ, abs.(collect(krange))]
end

let 
    data = res_area_prep(pop2)
    rescurve = figure()
    plot(data[1], data[2])
    plot(data[1], 80 .* tanh.(data[1] ./ 5)) #?????part of the proof - because res bound for real model is a subset of quad and we know quad has finite area
    xlabel("τ")
    ylabel("k")
    return rescurve
end

function res_area(res_curve_data)
    τ = reverse(res_curve_data[1])
    k = reverse(res_curve_data[2])
    area1=trapz(τ,80 .- k) #trap rule d(tau) with unequal intervals
    area2=trapz(k,τ) #trap rule d(kick) with equal intervals
    return (area1+area2)/2
end

res_area(res_area_prep(pop1))
res_area(res_area_prep(pop2))