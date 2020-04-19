
function instantiate_initial_conditions(model)# system::PSY.System)
        #TODO: SolvePowerFlow here for eg_d, eg_q and others if needed.
        _initial_guess = [
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
        # TODO: use Nlsolve here properly to get valid initial conditions
        initial_conditions = Array{Pair}(undef, length(_initial_guess))
        for (ix, val) in enumerate(_initial_guess)
            initial_conditions[ix] = MTK.states(model)[ix] => val
        end
    return initial_conditions
end
