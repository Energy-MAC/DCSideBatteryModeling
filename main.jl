using OrdinaryDiffEq #Gets the solvers
using Plots
using PowerSystems

include(joinpath(pwd(), "DCSideBatteryModeling", "DCSideBatteryModeling.jl"))
# Only need to run this line to re-generate the system data
#include(joinpath(pwd(), "data","make_data.jl"))
# Load Data with PF solution from file
omib_sys = System(joinpath(pwd(), "data", "OMIB_inverterDCside.json"))

_parameter_values = instantiate_parameters(omib_sys)
M = instantiate_model(omib_sys)
u0 = M(_parameter_values)

# The use of these methods causes StackOverflow
#jac = instantiate_jacobian(M)
#jac(M)

# WIP Jacobian Experiments. This is the function isntantiate evaluated in main to avoid
# world age errors
jac_exp = get_jacobian_expression()
_jac = eval(jac_exp)
jac_eval = (out, u0, params) -> _jac(out, u0, params)
param_eval = (out, params) -> _jac(out, M.u0, params)
n = length(M.u0)
J = zeros(n, n)
_parameter_values = [x.second for x in M.parameters]
param_eval(J, _parameter_values)
jac = ModelJacobian(jac_eval, J)
jac(M)

# Returns Generic ODE system and solves
ode_prob = instantiate_ode(omib_sys; tspan = (0.0, 5))
sol1 = solve(ode_prob, Rosenbrock23())
plot(sol1, vars = (0, 13), title = "DC Voltage Before Load Step")

#=
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
