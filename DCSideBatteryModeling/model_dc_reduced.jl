function reduced_dc_model(::Nothing)
    # Model Parameters
    params = MTK.@parameters begin
        t
        # AC side quantities
        ωb      # Base Frequency
        # Grid impadance
        lg      # Grid reactance
        rg      # Grid resistance
        #Reference set-point input
        pʳ      # Active Power Setpoint
        qʳ      # Reactive Power Setpoint
        vʳ      # Voltage Setpoint
        ωʳ      # Frequency Setpoint
        # Load at rated voltage
        vl      # Load rated voltage
        pl      # Load rated power
        # Filter parameters
        lf      # Filter reactance
        cf      # Filter capacitance
        rf      # Filter Resistance
        ωz      # Filtering frequency
        # Transformer Parameters
        rt      # Transformer resistance
        lt      # Transformer reactance
        # OuterControl Loops
        Dp      # Active Power Droop
        Dq      # Reactive Power Droop
        # SRF Current Control
        kip     # Current control propotional gain
        kii     # Current control integral gain
        kffi    # Current control differential gain
        # SRF Voltage Control
        kvp     # Voltage control propotional gain
        kvi     # Voltage control integral gain
        kffv    # Voltage control differential gain
        # Virtual Impedance
        rv
        lv
        # DC Source Parameters
        ldc     # Equivalent inductance (i.e. battery inductance and DC/DC converter inductance)
        req
        rb0     # Equivalent resistance
        lb1
        rl1
        lb2
        rl2
        cb1
        rc1
        cb2
        rc2
        vb      # Battery Voltage
        cdc     # Dc-side capacitance
        # DC/DC converter controller parameters
        vdcʳ    # DC Voltage reference
        kpvb    # DC/DC Voltage control integral gain
        kivb    # DC/DC Voltage control propotional gain
        kpib    # DC/DC Current control propotional gain
        kiib    # DC/DC Current control Integral gain
        kpred
        a1      # First coefficient of Pade approximation
        a2      # Second co-efficient of Pade approxmiation
        b1
        b2
        b3
        Ts      # DC/DC controller time delay
    end

    MTK.@derivatives d'~t

    # Definition of the states
    states = MTK.@variables begin
        vdc(t)  #DC Voltage
        ibat(t) #Battery Current
        il1(t)
        il2(t)
        vc1(t)
        vc2(t)
        η(t)    #Integrator term for outer DC/DC PI controller
        κ(t)    #Integrator term for inner DC/DC PI controller
        A(t)    # First term for Pade approx
        B(t)    # Second term for Pade approx.
    end

    # Definition of the variables for non-linear system. Requires https://github.com/SciML/ModelingToolkit.jl/issues/322 to eliminate
    variables = MTK.@variables begin
        vdc  #DC Voltage
        ibat #Battery Current
        il1
        il2
        vc1
        vc2
        η    #Integrator term for outer DC/DC PI controller
        κ    #Integrator term for inner DC/DC PI controller
        A    # First term for Pade approx
        B    # Second term for Pade approx
    end

    # Expressions
    #p_inv = v_md * is_d + v_mq * is_q # Active power drawn from inverter
    p_inv = pl#eg_d * is_d + eg_q * is_q # Active power drawn from inverter
    i_ref = kpvb * (vdcʳ - vdc) + kivb * η
    i_in = (vb * ibat - ibat^2 * rb0) / vdc
    d_dc =(a2/Ts) * A + kpib * (i_ref + p_inv / (vdc) - i_in) + kiib * κ 

    model_rhs = [
        ### DC-side equations
        #∂vdc/∂t
        (ωb/cdc) * (i_in - p_inv / (vdc))
        #∂ibat/∂t
        (ωb / ldc) * (vb - rl2*(ibat-il2) - rl1*(ibat-il1) - vc2 - vc1 - rb0 * ibat - (1 - d_dc) * vdc)
        #∂il1/∂t
        (ωb / lb1) * (rl1*(ibat-il1))
        #∂il2/∂t
        (ωb / lb2) * (rl2*(ibat-il2))
        #∂vc1/∂t
        (ωb / cb1) * ibat-vc1/rc1
        #∂vc2/∂t
        (ωb / cb2) * ibat-vc2/rc2
        #∂η/∂t
        vdcʳ - vdc # Integrator for DC/DC outer PI controller
        #∂κ/dt
        i_ref + p_inv / (vdc) - i_in # Integrator for DC/DC inner PI controller
        # ∂A/dt
        (a1/Ts)*A + (a2/Ts^2)*B + kpib * (i_ref + p_inv / (vdc) - i_in) + kiib * κ # First term in Pade approximation
        # ∂B/dt
        A # Second term in Pade approx.
    ]

    # Temporary until SteadyState problems are resolved.
    model_lhs = [
        d(vdc)
        d(ibat)
        d(il1)
        d(il2)
        d(vc1)
        d(vc2)
        d(η)
        d(κ)
        d(A)
        d(B)
    ]

    return model_lhs, model_rhs, states, variables, params
end

function get_reduced_dc_nonlinear_system()
    _, model_rhs, _, variables, params = reduced_dc_model(nothing)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    return MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
end
