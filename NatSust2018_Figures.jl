include("JuliaPackages.jl")
include("FlowKickTrajectories.jl")
include("ResilienceBoundary.jl")
include("Newton.jl")
include("Branch.jl")
include("ResilienceArea.jl")

#THIS SCRIPT CALCULATES FLOW-KICK TRAJECTORIES FOR A HARVESTED POPULATION
function allee!(du, u, p, t)
    du[1] = u[1] * (1.0 - (u[1]/100.0)) * ((u[1]/20.0) - 1.0)
end

let 
    data1 = flow_kick_trajectory(allee!,-12.0, 0.25, 100.0, 64)
    data2 = flow_kick_trajectory(allee!,-40.0, 10/12, 100.0, 15)
    flowkicktrajectory = figure()
    plot(data1.t, data1.u)
    plot(data2.t, data2.u)
    xlabel("Time (years)")
    ylabel("Population Size \n (kt biomass)")
    return flowkicktrajectory
end
#THIS SECTION PLOTS THE FLOW-KICK RESILIENCE BOUNDARY FOR ALLEE POPULATIONS 1 AND 2
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

#THIS SCRIPT PLOTS THE FLOW-KICK RESILIENCE BOUNDARY FOR A MODEL LAKE
#For each nutrient pulse (kick) size, the script approximates the minimum time it takes for phosphorus levels to recover by that magnitude. 

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


#This script traces out a branch of fixed points for a flow-kick map using Newton's method and parameter continuation

let 
    data1 = branch(0.0, 20.0, 0.001, 0.25, 100.0)
    data2 = branch(0.0, 20.0, 0.001, 0.25, 20.0)
    data3 = branch(0.0, 20.0, 0.001, 0.25, 0.0)
    test = figure()
    plot(data1[1], data1[2])
    plot(data2[1], data2[2])
    plot(data3[1], data3[2])
    return test
end


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
    datax00y00=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.0, 0.0)
    datax00y02=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.0, 0.2)
    datax00y04=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.0, 0.4)
    datax00y10=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.0, 1.0)
    datax10y05=flow_trajectory(Stommel!, StommelPar(), 200.0, 1.0, 0.5)
    datax01y00=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.1, 0.0)
    datax03y00=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.3, 0.0)
    datax07y00=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.7, 0.0)
    datax01y10=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.1, 1.0)
    datax03y10=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.3, 1.0)
    datax04y10=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.4, 1.0)
    datax045y10=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.45, 1.0)
    datax05y10=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.5, 1.0)
    datax07y10=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.7, 1.0)
    datax09y10=flow_trajectory(Stommel!, StommelPar(), 200.0, 0.9, 1.0)
    datasep28 = flow_trajectory(MinusStommel!, StommelPar(), 2.8, 0.35184, 0.7651-0.001)
    datasep17 = flow_trajectory(MinusStommel!, StommelPar(), 1.7, 0.35184, 0.7651+0.001)
    datastableman1 = flow_trajectory(Stommel!, StommelPar(), 200.0, 0.35184-0.001, 0.7651)
    datastableman2 = flow_trajectory(Stommel!, StommelPar(), 200.0, 0.35184+0.001, 0.7651)
    test = figure()
    plot(datax00y00[1,:], datax00y00[2,:], color="black")
    plot(datax00y02[1,:], datax00y02[2,:], color="black")
    plot(datax00y04[1,:], datax00y04[2,:], color="black")
    plot(datax00y10[1,:], datax00y10[2,:], color="black")
    plot(datax01y00[1,:], datax01y00[2,:], color="black")
    plot(datax03y00[1,:], datax03y00[2,:], color="black")
    plot(datax07y00[1,:], datax07y00[2,:], color="black")
    plot(datax01y10[1,:], datax01y10[2,:], color="black")
    plot(datax03y10[1,:], datax03y10[2,:], color="black")
    plot(datax04y10[1,:], datax04y10[2,:], color="black")
    plot(datax045y10[1,:], datax045y10[2,:], color="black")
    plot(datax05y10[1,:], datax05y10[2,:], color="black")
    plot(datax07y10[1,:], datax07y10[2,:], color="black")
    plot(datax09y10[1,:], datax09y10[2,:], color="black")
    plot(datasep28[1,:], datasep28[2,:], color="black", linestyle="dashed")
    plot(datasep17[1,:], datasep17[2,:], color="black", linestyle="dashed")
    plot(datastableman1[1,:], datastableman1[2,:], color="black")
    plot(datastableman2[1,:], datastableman2[2,:], color="black")
    plot(stablept1[1],stablept1[2], "bo")
    plot(stablept2[1],stablept2[2], "bo")
    plot(unstablept[1],unstablept[2], "bo")
    xlim(0,1)
    ylim(0,1)
    return test
end
    
#Exhibit outcome of short recovery period: stabilize in basin of attr. of a
#(Figure 5b)

let 
    stablept1 = [0.135, 0.48358]
    stablept2 =[0.43205, 0.82028]
    unstablept = [0.35184, 0.7651]
    datasep28 = flow_trajectory(MinusStommel!, StommelPar(), 2.8, 0.35184, 0.7651-0.001)
    datasep17 = flow_trajectory(MinusStommel!, StommelPar(), 1.7, 0.35184, 0.7651+0.001)
    data = flow_kick_trajectory(Stommel!, StommelPar(), 0.1, 0.1, [0.135, 0.48358], 20)
    flowkicktrajectory = figure()
    plot(data[1,:], data[2,:])
    plot(datasep28[1,:], datasep28[2,:], color="black", linestyle="dashed")
    plot(datasep17[1,:], datasep17[2,:], color="black", linestyle="dashed")
    plot(stablept1[1],stablept1[2], "bo")
    plot(stablept2[1],stablept2[2], "bo")
    plot(unstablept[1],unstablept[2], "bo")
    xlabel("Salinity Anomaly")
    ylabel("Temperature anomaly")
    xlim(0,1)
    ylim(0,1)
    return flowkicktrajectory
end


#Exhibit outcome of long recovery period: stabilize in basin of attr. of c
#(Figure 5c)
flow_kick_trajectory(model::Function, p, k, τ, x0, N)

let 
    stablept1 = [0.135, 0.48358]
    stablept2 =[0.43205, 0.82028]
    unstablept = [0.35184, 0.7651]
    datasep28 = flow_trajectory(MinusStommel!, StommelPar(), 2.8, 0.35184, 0.7651-0.001)
    datasep17 = flow_trajectory(MinusStommel!, StommelPar(), 1.7, 0.35184, 0.7651+0.001)
    data = flow_kick_trajectory(Stommel!, StommelPar(), 0.1, 1, [0.135, 0.48358], 10)
    flowkicktrajectory = figure()
    plot(data[1,:], data[2,:])
    plot(datasep28[1,:], datasep28[2,:], color="black", linestyle="dashed")
    plot(datasep17[1,:], datasep17[2,:], color="black", linestyle="dashed")
    plot(stablept1[1],stablept1[2], "bo")
    plot(stablept2[1],stablept2[2], "bo")
    plot(unstablept[1],unstablept[2], "bo")
    xlabel("Salinity Anomaly")
    ylabel("Temperature anomaly")
    xlim(0,1)
    ylim(0,1)
    return flowkicktrajectory
end

#Exhibit outcome of switching disturbance pattern to recovery time of 1 once stabilized at flow-kick equilibrium for recovery time 0.1 (Figure 5d)
let 
    stablept1 = [0.135, 0.48358]
    stablept2 =[0.43205, 0.82028]
    unstablept = [0.35184, 0.7651]
    datasep28 = flow_trajectory(MinusStommel!, StommelPar(), 2.8, 0.35184, 0.7651-0.001)
    datasep17 = flow_trajectory(MinusStommel!, StommelPar(), 1.7, 0.35184, 0.7651+0.001)
    data_1st = flow_kick_trajectory(Stommel!, StommelPar(), 0.1, 0.1, [0.135, 0.48358], 20)
    data_2nd = flow_kick_trajectory(Stommel!, StommelPar(), 0.1, 1, [data_1st[1,end], data_1st[2,end]], 20)
    flowkicktrajectory = figure()
    plot(data_1st[1,:], data_1st[2,:])
    plot(data_2nd[1,:], data_2nd[2,:])
    plot(datasep28[1,:], datasep28[2,:], color="black", linestyle="dashed")
    plot(datasep17[1,:], datasep17[2,:], color="black", linestyle="dashed")
    plot(stablept1[1],stablept1[2], "bo")
    plot(stablept2[1],stablept2[2], "bo")
    plot(unstablept[1],unstablept[2], "bo")
    xlabel("Salinity Anomaly")
    ylabel("Temperature anomaly")
    xlim(0,1)
    ylim(0,1)
    return flowkicktrajectory
end