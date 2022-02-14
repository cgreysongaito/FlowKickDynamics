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

function flow_trajectory(model::Function, p, τ, x0, y0)
    tspan = (0.0, τ)
    tvals = 0.0:0.1:τ
    prob = ODEProblem(model, [x0, y0], tspan, p)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    return solend = sol(tvals)
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