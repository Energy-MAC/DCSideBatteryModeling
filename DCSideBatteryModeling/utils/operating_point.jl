struct ModelOperatingPoint
    sys_func::Function
    u0::Vector{Float64}
    parameters::Vector
end

function instantiate_model(system::PSY.System; solve_powerflow = false)
    nl_sys = get_nonlinear_system()
    variable_count = length(states(nl_sys))
    nlsys_func = MTK.generate_function(nl_sys, expression = Val{false})[2]
    sys_f = (out, x, param) -> nlsys_func(out, x, param)

    if solve_powerflow
        PSY.solve_powerflow!(system, NLsolve.nlsolve)
    end

    inverter = collect(PSY.get_components(PSY.DynamicInverter, system))
    isempty(inverter) && @error("There are no inverters in the system")
    bus_voltage = PSY.get_voltage(PSY.get_bus(inverter[1]))
    bus_angle = PSY.get_angle(PSY.get_bus(inverter[1]))
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
        0.5,    #ibat
        0.0,    #η
        0.0,    #κ
        0.0,    #M
        0.0,    #L
        ]
    parameter_values = instantiate_parameters(system, nl_sys)
    M = ModelOperatingPoint(sys_f, initial_guess, parameter_values)
    M()
    return M
end

function (M::ModelOperatingPoint)(parameters::Vector)
    _parameter_values = [x.second for x in parameters]
    res = NLsolve.nlsolve((out, x) -> M.sys_func(out, x, _parameter_values), M.u0)
    !NLsolve.converged(res) && @error("NLsolve failed to converge")
    M.parameters .= parameters
    M.u0 .= res.zero
    return M.u0
end

(M::ModelOperatingPoint)() = M(M.parameters)
