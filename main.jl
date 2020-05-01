using OrdinaryDiffEq #Gets the solvers
using Plots
include(joinpath(pwd(), "DCSideBatteryModeling", "DCSideBatteryModeling.jl"))

# Returns Generic ODE system
model = get_model();
f = MTK.generate_function(model)[2]
ode_prob = instantiate_model(model, (0.0, 0.1))
sol1 = solve(ode_prob, TRBDF2()); #Use Solver for stiff problems
plot(sol1, vars = (0, 13), title = "DC Voltage Before Load Step")

#_parameter_values = instantiate_parameters(model) #, system)
#parameter_values = [x.second for x in _parameter_values]
#jac = get_jacobian_function();

# WIP Jacobian Experiments
_, model_rhs, _, variables, params = get_internal_model(nothing)
variable_count = length(variables)
_eqs = zeros(length(model_rhs)) .~ model_rhs
_nl_system = MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])

# This works
nlsys_jac = MTK.generate_jacobian(_nl_system)[2] # second is in-place
jac = eval(nlsys_jac)
param_eval = (out, params) -> jac(out, ode_prob.u0, params)
n= length(ode_prob.u0)
J = zeros(n, n)
param_eval(J, parameter_values)

# This causes StackOverflow works
jac = MTK.generate_jacobian(_nl_system, expression = Val{false})[2] # second is in-place
param_eval = (out, params) -> jac(out, ode_prob.u0, params)
n= length(ode_prob.u0)
J = zeros(n, n)
param_eval(J, parameter_values)
param_eval = (out, params) -> jac(out, ode_prob.u0, params)
n= length(ode_prob.u0)
J = zeros(n, n)
param_eval(J, parameter_values)

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
