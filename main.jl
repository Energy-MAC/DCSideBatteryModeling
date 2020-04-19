include(joinpath(pwd(), "DCSideBatteryModeling","DCSideBatteryModeling.jl"))

# Returns Generic ODE system
system = get_system()

#=
parameters = get_params(Ub,fb,ωb,Sb,Vb)

ic = get_IC()

tspan = (0.0,0.5)
prob = ODEProblem(ode_system!,ic.zero,tspan, parameters)
sol1 = solve(prob)

plot(sol1,vars=(0,13),title = "DC Voltage Before Load Step")

parameters.pl = 0.6;

tspan = (0.0,1)
prob = ODEProblem(ode_system!,sol1.u[end],tspan, parameters)
sol2 = solve(prob)

plot(sol2,vars=(0,13),title = "DC Voltage After Load Step")

function condition(u,t,integrator)
    t == 0.2
end

function affect!(integrator)
  parameters.pl = 0.6;
end
cb = DiscreteCallback(condition, affect!)

parameters = get_params(Ub,fb,ωb,Sb,Vb)
const tstop = [0.2]

parameters.pl

prob = ODEProblem(ode_system!,ic.zero,tspan, parameters)
sol3 = solve(prob,Tsit5(),callback = cb, tstops=tstop)

plot(sol1,vars=(0,13),title = "DC Voltage Before Load Step")

parameters.pl
=#
