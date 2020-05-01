function solve_steady_state(initial_guess, parameter_values)
    _, model_rhs, _, variables, params = get_internal_model(nothing)
    @assert length(initial_guess) == length(model_rhs) == length(variables)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    _nl_system = MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
    #_parameter_values = [x.second for x in parameter_values]
    nlsys_jac = MTK.generate_jacobian(_nl_system, expression = Val{false})[2] # second is in-place
    return nlsys_jac
end


    (out, x) -> nlsys_func(out, x, _parameter_values),
