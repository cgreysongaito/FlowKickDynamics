#THIS SCRIPT PLOTS THE FLOW-KICK RESILIENCE BOUNDARY FOR A MODEL LAKE
#For each nutrient pulse (kick) size, the script approximates the minimum time it takes for phosphorus levels to recover by that magnitude. 
include("JuliaPackages.jl")

#Parameters
@with_kw mutable struct LakePar
    l=25 #background input rate of P from the watershed
    s=0.5 #linear P loss rate
    r=50 #maximum P recycling rate from sediments
    q=8 #parametrizes sigmoid recycling curve shape
    m=100 #half saturation constant for recycling
end

function lake_model(x, p)
    @unpack l, s, r, q, m = p
    return l - s*x + (r*(x^q))/(m^q + x^q)
end

let
    data1 = res_bound(lake_model, LakePar(), 0.0:1:49.5, 50.0, 100.0, 0.5)
    test = figure()
    plot(data1[1],data1[2])
    xlabel("flow time (τ)")
    ylabel("kick size (κ)")
    xlim(0,15)
    ylim(0,50)
    return test
end