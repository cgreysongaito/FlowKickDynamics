include("JuliaPackages.jl")
include("FlowKickTrajectories.jl")
include("ResilienceBoundary.jl")
include("ResilienceArea.jl")
include("NewtonBranch.jl")


## THIS SECTION CALCULATES FLOW-KICK TRAJECTORIES FOR A HARVESTED POPULATION
#Parameters
@with_kw mutable struct AlleePar
    K=100.0  #carrying capacity
    A=20.0 #critical Allee threshold
    α=0.0 #first parameter in quadratic modifier
    β=0.0 #second parameter in quadratic modifier
    γ=1.0 #third parameter in quadratic modifier - set at 1.0 for when don't want quadratic modifier
end

function allee!(du, u, p, t,) #Used for producing time series (numerically solving the model)
    @unpack K, A, α, β, γ = p
    du[1] = (u[1] * (1.0 - (u[1]/K)) * ((u[1]/A) - 1.0)) * ((α * u[1]^2) - (β * u[1]) + γ)
    return
end

function allee(x, p) #used for finding the resilience boundary
    @unpack K, A, α, β, γ = p
    return (x * (1.0 - (x/K)) * ((x/A) - 1.0)) * ((α * x^2) - (β * x) + γ)
end

let 
    data1 = flow_kick_trajectory(allee!, AlleePar(), -12.0, 0.25, 64, [100.0])
    data2 = flow_kick_trajectory(allee!, AlleePar(), -40.0, 10/12, 15, [100.0])
    flowkicktrajectory = figure()
    plot(data1.t, data1.u, color="orange")
    plot(data2.t, data2.u, color="red")
    xlabel("Time (years)")
    ylabel("Population Size\n(kt biomass)")
    ylim(0, 110)
    return flowkicktrajectory
end #Figure 1b

## THIS SECTION PLOTS THE FLOW-KICK RESILIENCE BOUNDARY FOR ALLEE POPULATIONS 1 AND 2
let
    data1 = res_bound(allee, AlleePar(), -79.8:0.1:0.0, 100.0, 20.0, 0.1)
    data2 = res_bound(allee, AlleePar(α = 0.0002, β=0.024, γ=1.4), -79.8:0.1:0.0, 100.0, 20.0, 0.1)
    resbound = figure()
    plot(data1[1],data1[2], color = "blue")
    plot(data2[1],data2[2], color = "green", linestyle = "dashed")
    xlabel("Time interval between\nharvest (flow time)\n(years)")
    ylabel("Harvest size (kick)\n(kt biomass)")
    xlim(0,5)
    ylim(0,80)
    return resbound
end #Figure 1a

## THIS SECTION CALCULATES THE AREA OF THE NONRESILIENT SELECTION (AND PROVIDES GRAPH FOR PROOF THAT NONRESLIENT AREA FROM FULL MODEL IS SUBSET OF QUADRACTIC FUNCTION NONRESILIENT AREA)
#####WHY is res boundary calculated differently?
let 
    data = res_area_prep(allee, AlleePar(α = 0.0002, β=0.024, γ=1.4), -79.99:0.01:0.0, 100.0)
    rescurve = figure()
    plot(data[1], data[2])
    plot(data[1], 80 .* tanh.(data[1] ./ 5)) #?????part of the proof - because res bound for real model is a subset of quad and we know quad has finite area
    xlabel("τ")
    ylabel("k")
    return rescurve
end

res_area(res_area_prep(allee, AlleePar(), -79.99:0.01:0.0, 100.0), 80.0)
res_area(res_area_prep(allee, AlleePar(α = 0.0002, β=0.024, γ=1.4), -79.99:0.01:0.0, 100.0), 80.0)

## THIS SECTION PLOTS THE FLOW-KICK RESILIENCE BOUNDARY FOR THE LAKE ATER QUALITY MODEL
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
    resbound = figure()
    plot(data1[1],data1[2])
    xlabel("Recovery time between\npulses (months)")
    ylabel("Nutrient pulse\nmagnitude (hg)")
    xlim(0,15)
    ylim(0,50)
    return resbound
end #Figure 4b

## THIS SECTION TRACES OUT A BRANCH OF FIXED POINTS FOR A FLOW-KICK MAP USING NEWTON's METHOD AND PARAMETER CONTINUATION
function CoupledVar_lake!(du, u, p, t)
    @unpack l, s, r, q, m = p
    #gives ODE for fish population 1, along with varational equation to be solved simultaneously
    du[1] = l - s*u[1] + (r*(u[1]^q))/(m^q + u[1]^q)
    du[2] = (-s + ((m^q) * r * q * u[1]^(q-1))/(((m^q)+(u[1]^q))^2))*u[2]
end

let 
    data1 = branch(CoupledVar_lake!, LakePar(), 0.0:0.001:25.0, 2.0, 50.0)
    data2 = branch(CoupledVar_lake!, LakePar(), 0.0:0.001:25.0, 2.0, 100.0)
    data3 = branch(CoupledVar_lake!, LakePar(), 0.0:0.001:25.0, 2.0, 150.0)
    branchfig = figure()
    plot(data1[1], data1[2], color="black")
    plot(data1[1], data1[2]-data1[1], color="gray")
    plot(data2[1], data2[2], color="black", linestyle = "dashed")
    plot(data2[1], data2[2]-data2[1], color="gray", linestyle = "dashed")
    plot(data3[1], data3[2], color="black")
    plot(data3[1], data3[2]-data3[1], color="gray")
    return branchfig
end

function CoupledVar_allee!(du, u, p, t)
    @unpack K, A = p
    #gives ODE for fish population 1, along with varational equation to be solved simultaneously
    du[1] = u[1] * (1.0 - (u[1]/K)) * ((u[1]/A) - 1.0)
    du[2] = (-1+3*(u[1]/25) - 3*(u[1]^2/2000))*u[2]
end

let 
    data1 = branch(CoupledVar_allee!, AlleePar(), -20.0:0.001:0.0, 0.25, 100.0)
    data2 = branch(CoupledVar_allee!, AlleePar(), -20.0:0.001:0.0, 0.25, 20.0)
    data3 = branch(CoupledVar_allee!, AlleePar(), -20.0:0.001:0.0, 0.25, 0.0)
    branchfig = figure()
    plot(data1[1], data1[2], color="black")
    plot(data1[1], data1[2]-data1[1], color="gray")
    plot(data2[1], data2[2], color="black", linestyle = "dashed")
    plot(data2[1], data2[2]-data2[1], color="gray", linestyle = "dashed")
    plot(data3[1], data3[2], color="black")
    plot(data3[1], data3[2]-data3[1], color="gray")
    return branchfig
end #SI Figure 2b (but for the Allee model instead of the Lake model)

#THIS SCRIPT PLOTS THE PHASE PORTRAITS AND FLOW-KICK TRAJECTORIES FOR STOMMEL'S OCEAN BOX MODEL (FIGURE 5). 
#Although Stommel's model reduces to coordinates x and y corresponding to salinity and temperature in the low latitude box, the coordinates -x and -y were used in the paper to focus on the high latitude ocean box. to recreate Figure 5, simply rotate figures by 180 degrees and reverse the signs on the axes.

@with_kw mutable struct StommelPar
    λ=1/5  #carrying capacity
    R=2 #critical Allee threshold
    δ=1/6 #first parameter in quadratic modifier
end

function Stommel!(du, u, p, t,)
    @unpack λ, R, δ = p
    X, Y = u
    du[1]=δ*(1-X)-(X/λ)*abs(-Y+R*X)
    du[2]=1-Y-(Y/λ)*abs(-Y+R*X)
end

function MinusStommel!(du, u, p, t,)
    @unpack λ, R, δ = p
    X, Y = u
    du[1]=-(δ*(1-X)-(X/λ)*abs(-Y+R*X))
    du[2]=-(1-Y-(Y/λ)*abs(-Y+R*X))
end

#Create phase portrait (Figure 5a)
let 
    stablept1 = [0.135, 0.48358]
    stablept2 =[0.43205, 0.82028]
    unstablept = [0.35184, 0.7651]
    datax00y00 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.0, 0.0])
    datax00y02 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.0, 0.2])
    datax00y04 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.0, 0.4])
    datax00y10 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.0, 1.0])
    datax10y05 = flow_trajectory(Stommel!, StommelPar(), 200.0, [1.0, 0.5])
    datax01y00 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.1, 0.0])
    datax03y00 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.3, 0.0])
    datax07y00 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.7, 0.0])
    datax01y10 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.1, 1.0])
    datax03y10 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.3, 1.0])
    datax04y10 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.4, 1.0])
    datax045y10 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.45, 1.0])
    datax05y10 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.5, 1.0])
    datax07y10 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.7, 1.0])
    datax09y10 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.9, 1.0])
    datasep28 = flow_trajectory(MinusStommel!, StommelPar(), 2.8, [0.35184, 0.7651-0.001])
    datasep17 = flow_trajectory(MinusStommel!, StommelPar(), 1.7, [0.35184, 0.7651+0.001])
    datastableman1 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.35184-0.001, 0.7651])
    datastableman2 = flow_trajectory(Stommel!, StommelPar(), 200.0, [0.35184+0.001, 0.7651])
    phaseportrait = figure()
    plot(.-datax00y00[1,:], .-datax00y00[2,:], color="black")
    plot(.-datax00y02[1,:], .-datax00y02[2,:], color="black")
    plot(.-datax00y04[1,:], .-datax00y04[2,:], color="black")
    plot(.-datax00y10[1,:], .-datax00y10[2,:], color="black")
    plot(.-datax01y00[1,:], .-datax01y00[2,:], color="black")
    plot(.-datax03y00[1,:], .-datax03y00[2,:], color="black")
    plot(.-datax07y00[1,:], .-datax07y00[2,:], color="black")
    plot(.-datax01y10[1,:], .-datax01y10[2,:], color="black")
    plot(.-datax03y10[1,:], .-datax03y10[2,:], color="black")
    plot(.-datax04y10[1,:], .-datax04y10[2,:], color="black")
    plot(.-datax045y10[1,:], .-datax045y10[2,:], color="black")
    plot(.-datax05y10[1,:], .-datax05y10[2,:], color="black")
    plot(.-datax07y10[1,:], .-datax07y10[2,:], color="black")
    plot(.-datax09y10[1,:], .-datax09y10[2,:], color="black")
    plot(.-datasep28[1,:], .-datasep28[2,:], color="black", linestyle="dashed")
    plot(.-datasep17[1,:], .-datasep17[2,:], color="black", linestyle="dashed")
    plot(.-datastableman1[1,:], .-datastableman1[2,:], color="black")
    plot(.-datastableman2[1,:], .-datastableman2[2,:], color="black")
    plot(-stablept1[1],-stablept1[2], "ko")
    plot(-stablept2[1],-stablept2[2], "ko")
    scatter(-unstablept[1],-unstablept[2], c="w", alpha=1, edgecolors="k")
    xlim(-1, 0)
    ylim(-1, 0)
    xlabel("Salinity Anomaly (non-dimensionalized)")
    ylabel("Temperature anomaly (non-dimensionalized)")
    return phaseportrait
end
    
#Exhibit outcome of short recovery period: stabilize in basin of attr. of a
#(Figure 5b)

let 
    stablept1 = [0.135, 0.48358]
    stablept2 =[0.43205, 0.82028]
    unstablept = [0.35184, 0.7651]
    datasep28 = flow_trajectory(MinusStommel!, StommelPar(), 2.8, [0.35184, 0.7651-0.001])
    datasep17 = flow_trajectory(MinusStommel!, StommelPar(), 1.7, [0.35184, 0.7651+0.001])
    data = flow_kick_trajectory(Stommel!, StommelPar(), 0.1, 0.1, 20, [0.135, 0.48358])
    flowkicktrajectory = figure()
    plot(.-data[1,:], .-data[2,:])
    plot(.-datasep28[1,:], .-datasep28[2,:], color="black", linestyle="dashed")
    plot(.-datasep17[1,:], .-datasep17[2,:], color="black", linestyle="dashed")
    plot(-stablept1[1],-stablept1[2], "ko")
    plot(-stablept2[1],-stablept2[2], "ko")
    scatter(-unstablept[1],-unstablept[2], c="w", alpha=1, edgecolors="k")
    xlabel("Salinity Anomaly (non-dimensionalized)")
    ylabel("Temperature anomaly (non-dimensionalized)")
    xlim(-1, 0)
    ylim(-1, 0)
    return flowkicktrajectory
end


#Exhibit outcome of long recovery period: stabilize in basin of attr. of c
#(Figure 5c)
let 
    stablept1 = [0.135, 0.48358]
    stablept2 =[0.43205, 0.82028]
    unstablept = [0.35184, 0.7651]
    datasep28 = flow_trajectory(MinusStommel!, StommelPar(), 2.8, [0.35184, 0.7651-0.001])
    datasep17 = flow_trajectory(MinusStommel!, StommelPar(), 1.7, [0.35184, 0.7651+0.001])
    data = flow_kick_trajectory(Stommel!, StommelPar(), 0.1, 1, 10, [0.135, 0.48358])
    flowkicktrajectory = figure()
    plot(.-data[1,:], .-data[2,:], color="red")
    plot(.-datasep28[1,:], .-datasep28[2,:], color="black", linestyle="dashed")
    plot(.-datasep17[1,:], .-datasep17[2,:], color="black", linestyle="dashed")
    plot(-stablept1[1],-stablept1[2], "ko")
    plot(-stablept2[1],-stablept2[2], "ko")
    scatter(-unstablept[1],-unstablept[2], c="w", alpha=1, edgecolors="k")
    xlabel("Salinity Anomaly (non-dimensionalized)")
    ylabel("Temperature anomaly (non-dimensionalized)")
    xlim(-1, 0)
    ylim(-1, 0)
    return flowkicktrajectory
end

#Exhibit outcome of switching disturbance pattern to recovery time of 1 once stabilized at flow-kick equilibrium for recovery time 0.1 (Figure 5d)
let 
    stablept1 = [0.135, 0.48358]
    stablept2 =[0.43205, 0.82028]
    unstablept = [0.35184, 0.7651]
    datasep28 = flow_trajectory(MinusStommel!, StommelPar(), 2.8, [0.35184, 0.7651-0.001])
    datasep17 = flow_trajectory(MinusStommel!, StommelPar(), 1.7, [0.35184, 0.7651+0.001])
    data_1st = flow_kick_trajectory(Stommel!, StommelPar(), 0.1, 0.1, 20, [0.135, 0.48358])
    data_2nd = flow_kick_trajectory(Stommel!, StommelPar(), 0.1, 1, 20, [data_1st[1,end], data_1st[2,end]])
    flowkicktrajectory = figure()
    plot(.-data_1st[1,:], .-data_1st[2,:], color="blue")
    plot(.-data_2nd[1,:], .-data_2nd[2,:], color="red")
    plot(.-datasep28[1,:], .-datasep28[2,:], color="black", linestyle="dashed")
    plot(.-datasep17[1,:], .-datasep17[2,:], color="black", linestyle="dashed")
    plot(-stablept1[1],-stablept1[2], "ko")
    plot(-stablept2[1],-stablept2[2], "ko")
    scatter(-unstablept[1],-unstablept[2], c="w", alpha=1, edgecolors="k")
    xlabel("Salinity Anomaly (non-dimensionalized)")
    ylabel("Temperature anomaly (non-dimensionalized)")
    xlim(-1, 0)
    ylim(-1, 0)
    return flowkicktrajectory
end