function instantiate_4th_order_dae(parameters, M::ModelOperatingPoint; tspan::Tuple,)
    model = get_4th_order_dae_system()
    initial_conditions = Array{Pair}(undef, length(M.u0))
    parameter_mapping = Array{Pair}(undef, length(M.parameters))
    for (ix, val) in enumerate(M.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    for (ix, val) in enumerate(parameters)
        parameter_mapping[ix] = MTK.parameters(model)[ix] => val
    end
    Mass = calculate_massmatrix(model)
    f = ODEFunction(model)

    return DiffEqBase.ODEProblem(
        f,
        [M.u0; 0.5; 0.5; 0.0],
        tspan,
        parameters,
        jac = false,
    )
#     return DiffEqBase.ODEProblem(
#         model,
#         initial_conditions,
#         tspan,
#         parameter_mapping,
#         jac = false,
#     )
end

function instantiate_2nd_order_dae(parameters, M::ModelOperatingPoint; tspan::Tuple,)
    model = get_2nd_order_dae_system()
    initial_conditions = Array{Pair}(undef, length(M.u0))
    parameter_mapping = Array{Pair}(undef, length(M.parameters))
    for (ix, val) in enumerate(M.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    for (ix, val) in enumerate(parameters)
        parameter_mapping[ix] = MTK.parameters(model)[ix] => val
    end
    Mass = calculate_massmatrix(model)
    f = ODEFunction(model)

    return DiffEqBase.ODEProblem(
        f,
        [M.u0; 0.5; 0.5; 0.0],
        tspan,
        parameters,
        jac = false,
    )
#     return DiffEqBase.ODEProblem(
#         model,
#         initial_conditions,
#         tspan,
#         parameter_mapping,
#         jac = false,
#     )
end

function instantiate_0th_order_dae(parameters, M::ModelOperatingPoint; tspan::Tuple,)
    model = get_0th_order_dae_system()
    initial_conditions = Array{Pair}(undef, length(M.u0))
    parameter_mapping = Array{Pair}(undef, length(M.parameters))
    for (ix, val) in enumerate(M.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    for (ix, val) in enumerate(parameters)
        parameter_mapping[ix] = MTK.parameters(model)[ix] => val
    end
    Mass = calculate_massmatrix(model)
    f = ODEFunction(model)

    return DiffEqBase.ODEProblem(
        f,
        [M.u0; 0.5; 0.5; 0.0],
        tspan,
        parameters,
        jac = false,
    )
#     return DiffEqBase.ODEProblem(
#         model,
#         initial_conditions,
#         tspan,
#         parameter_mapping,
#         jac = false,
#     )
end

function instantiate_4th_order_ode(system::PSY.System; tspan::Tuple, kwargs...)
    model = ode_model_4th_order()
    steady_state = instantiate_model(system; solve_powerflow = true)
    initial_conditions = Array{Pair}(undef, length(steady_state.u0))
    for (ix, val) in enumerate(steady_state.u0)
        initial_conditions[ix] = MTK.states(model)[ix] => val
    end
    
    return DiffEqBase.ODEProblem(
        model,
        initial_conditions,
        tspan,
        steady_state.parameters,
        jac = false,
    )
end
