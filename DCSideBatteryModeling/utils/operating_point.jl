struct ModelOperatingPoint
    sys_func::Function
    u0::Vector{Float64}
    parameters::Vector{Float64}
end


function instantiate_4th_order_model(system::PSY.System; solve_powerflow = false)
    nl_sys = get_4th_order_nonlinear_system()
    variable_count = length(states(nl_sys))
    nlsys_func = MTK.generate_function(nl_sys, expression = Val{false})[2]
    sys_f = (out, x, param) -> nlsys_func(out, x, param)

    if solve_powerflow
        PSY.solve_powerflow!(system, NLsolve.nlsolve)
    end

    #inverter = collect(PSY.get_components(PSY.DynamicInverter, system))
    #isempty(inverter) && @error("There are no inverters in the system")
    #bus_voltage = PSY.get_voltage(PSY.get_bus(inverter[1]))
    #bus_angle = PSY.get_angle(PSY.get_bus(inverter[1]))
    # TODO: instantiate initial guess from PSY power flow solution
    initial_guess = [
        1.0,    #eg_d
        0.0,    #eg_q
        0.5,    #is_d
        0,      #is_q
        0.5,    #ig_d
        0,      #ig_q
        0.5,
        0.0,
        0.0,    #ξ_d
        0.0,    #ξ_q
        0.0,    #γ_d
        0.0,    #γ_q
        1.0,    #vdc
        1.0,    #ibat
        1.0,    #il1
        1.0,    #il2
        0.0,
        0.0,
        0.0,    #η
        0.0,    #κ
        0.0,    #A
        0.0,    #B
        0.0,    #C
        0.0,    #D
        0.0,    #E
    ]
    parameter_mapping = instantiate_parameters(system, nl_sys)
    parameter_values = [x.second for x in parameter_mapping]
    M = ModelOperatingPoint(sys_f, initial_guess, parameter_values)
    M()
    return M
end

function instantiate_2nd_order_model(system::PSY.System; solve_powerflow = false)
    nl_sys = get_2nd_order_nonlinear_system()
    variable_count = length(states(nl_sys))
    nlsys_func = MTK.generate_function(nl_sys, expression = Val{false})[2]
    sys_f = (out, x, param) -> nlsys_func(out, x, param)

    if solve_powerflow
        PSY.solve_powerflow!(system, NLsolve.nlsolve)
    end

    #inverter = collect(PSY.get_components(PSY.DynamicInverter, system))
    #isempty(inverter) && @error("There are no inverters in the system")
    #bus_voltage = PSY.get_voltage(PSY.get_bus(inverter[1]))
    #bus_angle = PSY.get_angle(PSY.get_bus(inverter[1]))
    # TODO: instantiate initial guess from PSY power flow solution
    initial_guess = [
        1.0,    #eg_d
        0.0,    #eg_q
        0.5,    #is_d
        0,      #is_q
        0.5,    #ig_d
        0,      #ig_q
        0.5,
        0.0,
        0.0,    #ξ_d
        0.0,    #ξ_q
        0.0,    #γ_d
        0.0,    #γ_q
        1.0,    #vdc
        1.0,    #ibat
        0.01,
        0.01,
        0.0,    #η
        0.0,    #κ
        0.0,    #A
        0.0,    #B
        0.0,    #C
        0.0,    #D
        0.0,    #E
    ]
    parameter_mapping = instantiate_parameters(system, nl_sys)
    parameter_values = [x.second for x in parameter_mapping]
    M = ModelOperatingPoint(sys_f, initial_guess, parameter_values)
    M()
    return M
end

function instantiate_0th_order_model(system::PSY.System; solve_powerflow = false)
    nl_sys = get_0th_order_nonlinear_system()
    variable_count = length(states(nl_sys))
    nlsys_func = MTK.generate_function(nl_sys, expression = Val{false})[2]
    sys_f = (out, x, param) -> nlsys_func(out, x, param)

    if solve_powerflow
        PSY.solve_powerflow!(system, NLsolve.nlsolve)
    end

    #inverter = collect(PSY.get_components(PSY.DynamicInverter, system))
    #isempty(inverter) && @error("There are no inverters in the system")
    #bus_voltage = PSY.get_voltage(PSY.get_bus(inverter[1]))
    #bus_angle = PSY.get_angle(PSY.get_bus(inverter[1]))
    # TODO: instantiate initial guess from PSY power flow solution
    initial_guess = [
        1.0,    #eg_d
        0.0,    #eg_q
        0.5,    #is_d
        0,      #is_q
        0.5,    #ig_d
        0,      #ig_q
        0.5,
        0.0,
        0.0,    #ξ_d
        0.0,    #ξ_q
        0.0,    #γ_d
        0.0,    #γ_q
        1.0,    #vdc
        1.0,    #ibat
        0.0,    #η
        0.0,    #κ
        0.0,    #A
        0.0,    #B
        0.0,    #C
        0.0,    #D
        0.0,    #E
    ]
    parameter_mapping = instantiate_parameters(system, nl_sys)
    parameter_values = [x.second for x in parameter_mapping]
    M = ModelOperatingPoint(sys_f, initial_guess, parameter_values)
    M()
    return M
end

function (M::ModelOperatingPoint)(parameters::Vector{Float64})
    res = NLsolve.nlsolve((out, x) -> M.sys_func(out, x, parameters), M.u0, ftol=1e-8)
    !NLsolve.converged(res) && @error("NLsolve failed to converge")
    M.parameters .= parameters
    M.u0 .= res.zero
    return M.u0
end

function (M::ModelOperatingPoint)(
    parameters::Array{Pair{Variable{ModelingToolkit.Parameter{Number}}, Float64}, 1},
)
    parameter_values = [x.second for x in parameters]
    return M(parameter_values)
end

(M::ModelOperatingPoint)() = M(M.parameters)


# function (M::ModelOperatingPoint)(parameters::Vector)
#     _parameter_values = [x.second for x in parameters]
#     res = NLsolve.nlsolve((out, x) -> M.sys_func(out, x, _parameter_values), M.u0)
#     !NLsolve.converged(res) && @error("NLsolve failed to converge")
#     M.parameters .= parameters
#     M.u0 .= res.zero
#     return M.u0
# end

# (M::ModelOperatingPoint)() = M(M.parameters)



function solve_steady_state(initial_guess, parameter_values)
    _, model_rhs, _, variables, params = get_reduced_dc_model(nothing)
    @assert length(initial_guess) == length(model_rhs) == length(variables)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    _nl_system = MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
    nlsys_func = MTK.generate_function(_nl_system, expression = Val{false})[2]
    _parameter_values = [x.second for x in parameter_values]
    nlsys_jac = MTK.generate_jacobian(_nl_system, expression = Val{false})[2] # second is in-place
    sol = NLsolve.nlsolve(
        (out, x) -> nlsys_func(out, x, _parameter_values),
        #(out, x) -> nlsys_jac(out, x, _parameter_values),
        initial_guess,
    )
    println(sol)
    return sol.zero
end
function instantiate_initial_conditions(model, parameter_values)# system::PSY.System)
    #TODO: SolvePowerFlow here for eg_d, eg_q and others if needed.
    _initial_guess = [
        1.0,    #vdc
        0.5,    #ibat
        0.0,    #η
        0.0,    #κ
        0.0,    #M
        0.0,    #L
    ]
    _initial_conditions = solve_steady_state(_initial_guess, parameter_values)
    initial_conditions = Array{Pair}(undef, length(_initial_conditions))
    for (ix, val) in enumerate(_initial_conditions)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    return initial_conditions
end

function get_model()
    model_lhs, model_rhs, states, _, params = get_reduced_dc_model(nothing)
    t = params[1]
    return MTK.ODESystem(model_lhs .~ model_rhs, t, [states...], [params...][2:end])
end

function instantiate_model(
    model,
    tspan::Tuple,
    #system::PSY.System,
)
    parameter_values = instantiate_parameters(model) #, system)
    initial_conditions = instantiate_initial_conditions(model, parameter_values) #, system)
    #return DiffEqBase.ODEProblem(
    #    model,
    #    initial_conditions,
    #    tspan,
    #    parameter_values,
    #    jac = true,
    #)
    return initial_conditions
end