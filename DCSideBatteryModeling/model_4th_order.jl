function ode_model_4th_order(::Nothing)
    # Model Parameters
    params = MTK.@parameters begin
        t
        # Base quantities
        ωb      # Base Frequency
        Cb_dc   # DC-side base capacitance
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
        rv      # Virtual Resistance
        lv      # Vistual Inductance 
        # DC Source Parameters
        ldc     # Equivalent inductance (i.e. battery inductance and DC/DC converter inductance)
        req     # Equivalent resistance for 0th-order model
        rb0     # Full order resistor element 
        lb1     # Inductance of first RL branch
        rl1     # Resistance of first RL branch
        lb2     # Inductance of second RL branch  
        rl2     # Resistance of second RL branch
        cb1     # Capacitance of first RC branch
        rc1     # Inductance of first RC branch
        cb2     # Inductance of first RC branch
        rc2     # Inductance of first RC branch
        vb      # Battery Voltage
        cdc     # Dc-side capacitance
        # DC/DC converter controller parameters
        vdcʳ    # DC Voltage reference
        kpvb    # DC/DC Voltage control integral gain
        kivb    # DC/DC Voltage control propotional gain
        kpib    # DC/DC Current control propotional gain
        kiib    # DC/DC Current control Integral gain
        kpred   # DC/DC Current control predictive gain
        a1      # First coefficient of 2nd-order Pade approximation
        a2      # Second co-efficient of 2nd-order Pade approximation
        b1      # First coefficient of 3rd-order Pade approximation
        b2      # Second coefficient of 3rd-order Pade approximation
        b3      # Third coefficient of 3rd-order Pade approximation
        Ts      # DC/DC controller time delay
    end

    MTK.@derivatives d'~t

    # Definition of the states
    states = MTK.@variables begin
        eg_d(t) #d-axis capacitor filter voltage
        eg_q(t) #q-axis capacitor filter voltage
        is_d(t) #d-axis current flowing into filter
        is_q(t) #q-axis current flowing into filter
        ig_d(t) #d-axis current flowing into grid
        ig_q(t) #q-axis current flowing into grid
        pf(t)   #Filtered active power measurement
        qf(t)   #Filtered reactive power measurement
        ξ_d(t)  #d-axis integrator term for outer AC/DC PI controller
        ξ_q(t)  #q-axis integrator term for outer AC/DC PI controller
        γ_d(t)  #d-axis integrator term for inner AC/DC PI controller
        γ_q(t)  #d-axis integrator term for inner AC/DC PI controller
        vdc(t)  #DC Voltage
        ibat(t) #Battery Current
        il1(t)  #Current throught inductor in first RL branch 
        il2(t)  #Current throught inductor in second RL branch 
        vc1(t)  #Voltage across capacitor in first RC branch 
        vc2(t)  #Voltage across capacitor in second RC branch 
        η(t)    #Integrator term for outer DC/DC PI controller
        κ(t)    #Integrator term for inner DC/DC PI controller
        A(t)    #First term for 2nd-order Pade approx
        B(t)    #Second term for 2nd-order Pade approx.
        C(t)    #First term for 3rd-order Pade approx
        D(t)    #Second term for 3rd-order Pade approx.
        E(t)    #Third term for 3rd-order Pade approx
    end

    # Definition of the variables for non-linear system. Requires https://github.com/SciML/ModelingToolkit.jl/issues/322 to eliminate
    variables = MTK.@variables begin
        eg_d #d-axis capacitor filter voltage
        eg_q #q-axis capacitor filter voltage
        is_d #d-axis current flowing into filter
        is_q #q-axis current flowing into filter
        ig_d #d-axis current flowing into grid
        ig_q #q-axis current flowing into grid
        pf   #Filtered active power measurement
        qf   #Filtered reactive power measurement
        ξ_d  #d-axis integrator term for outer AC/DC PI controller
        ξ_q  #q-axis integrator term for outer AC/DC PI controller
        γ_d  #d-axis integrator term for inner AC/DC PI controller
        γ_q  #d-axis integrator term for inner AC/DC PI controller
        vdc  #DC Voltage
        ibat #Battery Current
        il1  #Current throught inductor in first RL branch 
        il2  #Current throught inductor in second RL branch 
        vc1  #Voltage across capacitor in first RC branch 
        vc2  #Voltage across capacitor in second RC branch 
        η    #Integrator term for outer DC/DC PI controller
        κ    #Integrator term for inner DC/DC PI controller
        A    # First term for Pade approx
        B    # Second term for Pade approx
        C    # First term for Pade approx
        D    # Second term for Pade approx
        E    # Second term for Pade approx
    end

    # Expressions
    pm = eg_d * ig_d + eg_q * ig_q  # AC Active Power Calculation
    qm = -eg_d * ig_q + eg_q * ig_d # AC Reactive Power Calculation
    ω_a = ωʳ + Dp * (pʳ - pf)  # Active Power Drop
    v_hat = vʳ + Dq * (qʳ - qf) # Reactive Power Drop
    v_iref_d = v_hat - rv * ig_d + ω_a * lv * ig_q # d-axis virtual impedance equation
    v_iref_q = -rv * ig_q - ω_a * lv * ig_d # q-axis virtual impedance equation
    i_hat_d = kvp * (v_iref_d - eg_d) + kvi * ξ_d - ω_a * cf * eg_q # Inner voltage controller d PI
    i_hat_q = kvp * (v_iref_q - eg_q) + kvi * ξ_q + ω_a * cf * eg_d # Inner voltage controller q PI
    v_md = kip * (i_hat_d - is_d) + kii * γ_d - ω_a * lf * is_q # Inner current controller d PI
    v_mq = kip * (i_hat_q - is_q) + kii * γ_q + ω_a * lf * is_d # Inner current controller q PI
    #p_inv = v_md * is_d + v_mq * is_q # Active power drawn from inverter
    p_inv = eg_d * is_d + eg_q * is_q # Active power drawn from inverter
    q_inv = -eg_d * is_q + eg_q * is_d # Reactive power drawn from inverter
    v_gd = (vl^2 / pl) * ig_d # d-qaxis Load voltage
    v_gq = (vl^2 / pl) * ig_q # q-axis load voltage
    i_ref = kpvb * (vdcʳ - vdc) + kivb * η # Current referecne for inner DC/DC PI controller 
    i_in = (vb * ibat - ibat^2 * rb0) / vdc # Current flwoing into DC-link capacitance 
    d_dc = (1 / 2) * ((a2 / Ts) * A) + (1 / 2) * ((-2 * b1 / Ts) * C + (-2 * b3 / Ts^3) * E) # DC/DC converter duty cycle
    is_d_dot = (ωb / lf) * (v_md - eg_d) - (rf * ωb / lf) * is_d + ωb * ω_a * is_q # derivative of is_d
    is_q_dot = (ωb / lf) * (v_mq - eg_q) - (rf * ωb / lf) * is_q - ωb * ω_a * is_d # derivative of is_q
    i_pred = (Ts / vdc) * (v_md * is_d_dot + v_mq * is_q_dot) # Estiamtion of change in current flowing out of DC-link capacitor over DC/DC converter switching period

    model_rhs = [
        ### Grid forming equations
        #𝜕eg_d/𝜕t
        (ωb / cf) * (is_d - ig_d) + ω_a * ωb * eg_q
        #𝜕eg_q/𝜕t
        (ωb / cf) * (is_q - ig_q) - ω_a * ωb * eg_d
        #𝜕is_d/𝜕t
        (ωb / lf) * (v_md - eg_d) - (rf * ωb / lf) * is_d + ωb * ω_a * is_q
        #𝜕is_q/𝜕t
        (ωb / lf) * (v_mq - eg_q) - (rf * ωb / lf) * is_q - ωb * ω_a * is_d
        #𝜕ig_d/𝜕t
        (ωb / lt) * (eg_d - v_gd) - (rt * ωb / lt) * ig_d + ωb * ω_a * ig_q
        #𝜕ig_q/𝜕t
        (ωb / lt) * (eg_q - v_gq) - (rt * ωb / lt) * ig_q - ωb * ω_a * ig_d
        #𝜕pf/𝜕t
        ωz * (pm - pf)
        #𝜕qf/𝜕t
        ωz * (qm - qf)
        #𝜕ξ_d/𝜕t
        v_iref_d - eg_d
        #𝜕ξ_q/𝜕t
        v_iref_q - eg_q
        #𝜕γ_d/𝜕t
        i_hat_d - is_d
        #𝜕γ_q/𝜕t
        i_hat_q - is_q
        ### DC-side equations
        #∂vdc/∂t
        (ωb / cdc) * (i_in - p_inv / (vdc))
        #∂ibat/∂t
        (ωb / ldc) * (
            vb - rl2 * (ibat - il2) - rl1 * (ibat - il1) - vc2 - vc1 - rb0 * ibat -
            (1 - d_dc) * vdc
        )
        #∂il1/∂t
        (ωb / lb1) * (rl1 * (ibat - il1))
        #∂il2/∂t
        (ωb / lb2) * (rl2 * (ibat - il2))
        #∂vc1/∂t
        (ωb / cb1) * (ibat - vc1 / rc1)
        #∂vc2/∂t
        (ωb / cb2) * (ibat - vc2 / rc2)
        #∂η/∂t
        vdcʳ - vdc # Integrator for DC/DC outer PI controller
        #∂κ/dt
        i_ref + p_inv / (vdc) - i_in # Integrator for DC/DC inner PI controller
        # ∂A/dt
        (a1 / Ts) * A +
        (a2 / Ts^2) * B +
        kpib * (i_ref + p_inv / (vdc) - i_in) +
        kiib * κ +
        kpred * i_pred  # First term in Pade approximation
        # ∂B/dt
        A # Second term in Pade approx.
        # ∂C/dt
        (b1 / Ts) * C +
        (b2 / Ts^2) * D +
        (b3 / Ts^3) * E +
        kpib * (i_ref + p_inv / (vdc) - i_in) +
        kiib * κ +
        kpred * i_pred  # First term in Pade approximation
        # ∂D/dt
        C
        # ∂E/dt
        D
    ]

    # Temporary until SteadyState problems are resolved.
    model_lhs = [
        d(eg_d)
        d(eg_q)
        d(is_d)
        d(is_q)
        d(ig_d)
        d(ig_q)
        d(pf)
        d(qf)
        d(ξ_d)
        d(ξ_q)
        d(γ_d)
        d(γ_q)
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
        d(C)
        d(D)
        d(E)
    ]

    return model_lhs, model_rhs, states, variables, params
end

function dae_model_4th_order(::Nothing)
    # Model Parameters
    params = MTK.@parameters begin
        t
        # Base quantities
        ωb      # Base Frequency
        Cb_dc   # DC-side base capacitance
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
        rv      # Virtual Resistance
        lv      # Vistual Inductance 
        # DC Source Parameters
        ldc     # Equivalent inductance (i.e. battery inductance and DC/DC converter inductance)
        req     # Equivalent resistance for 0th-order model
        rb0     # Full order resistor element 
        lb1     # Inductance of first RL branch
        rl1     # Resistance of first RL branch
        lb2     # Inductance of second RL branch  
        rl2     # Resistance of second RL branch
        cb1     # Capacitance of first RC branch
        rc1     # Inductance of first RC branch
        cb2     # Inductance of first RC branch
        rc2     # Inductance of first RC branch
        vb      # Battery Voltage
        cdc     # Dc-side capacitance
        # DC/DC converter controller parameters
        vdcʳ    # DC Voltage reference
        kpvb    # DC/DC Voltage control integral gain
        kivb    # DC/DC Voltage control propotional gain
        kpib    # DC/DC Current control propotional gain
        kiib    # DC/DC Current control Integral gain
        kpred   # DC/DC Current control predictive gain
        a1      # First coefficient of 2nd-order Pade approximation
        a2      # Second co-efficient of 2nd-order Pade approximation
        b1      # First coefficient of 3rd-order Pade approximation
        b2      # Second coefficient of 3rd-order Pade approximation
        b3      # Third coefficient of 3rd-order Pade approximation
        Ts      # DC/DC controller time delay
    end

    MTK.@derivatives d'~t

    # Definition of the states
    states = MTK.@variables begin
        eg_d(t) #d-axis capacitor filter voltage
        eg_q(t) #q-axis capacitor filter voltage
        is_d(t) #d-axis current flowing into filter
        is_q(t) #q-axis current flowing into filter
        ig_d(t) #d-axis current flowing into grid
        ig_q(t) #q-axis current flowing into grid
        pf(t)   #Filtered active power measurement
        qf(t)   #Filtered reactive power measurement
        ξ_d(t)  #d-axis integrator term for outer AC/DC PI controller
        ξ_q(t)  #q-axis integrator term for outer AC/DC PI controller
        γ_d(t)  #d-axis integrator term for inner AC/DC PI controller
        γ_q(t)  #d-axis integrator term for inner AC/DC PI controller
        vdc(t)  #DC Voltage
        ibat(t) #Battery Current
        il1(t)  #Current throught inductor in first RL branch 
        il2(t)  #Current throught inductor in second RL branch 
        vc1(t)  #Voltage across capacitor in first RC branch 
        vc2(t)  #Voltage across capacitor in second RC branch 
        η(t)    #Integrator term for outer DC/DC PI controller
        κ(t)    #Integrator term for inner DC/DC PI controller
        A(t)    #First term for 2nd-order Pade approx
        B(t)    #Second term for 2nd-order Pade approx.
        C(t)    #First term for 3rd-order Pade approx
        D(t)    #Second term for 3rd-order Pade approx.
        E(t)    #Third term for 3rd-order Pade approx
        d_dc(t) #DC/DC duty cycle
        i_in(t) #Current flowing into DC-link capacitance
        v_conv(t) #Magnitude of VSC output voltage
    end

    # Definition of the variables for non-linear system. Requires https://github.com/SciML/ModelingToolkit.jl/issues/322 to eliminate
    variables = MTK.@variables begin
        eg_d   #d-axis capacitor filter voltage
        eg_q   #q-axis capacitor filter voltage
        is_d   #d-axis current flowing into filter
        is_q   #q-axis current flowing into filter
        ig_d   #d-axis current flowing into grid
        ig_q   #q-axis current flowing into grid
        pf     #Filtered active power measurement
        qf     #Filtered reactive power measurement
        ξ_d    #d-axis integrator term for outer AC/DC PI controller
        ξ_q    #q-axis integrator term for outer AC/DC PI controller
        γ_d    #d-axis integrator term for inner AC/DC PI controller
        γ_q    #d-axis integrator term for inner AC/DC PI controller
        vdc    #DC Voltage
        ibat   #Battery Current
        il1    #Current throught inductor in first RL branch 
        il2    #Current throught inductor in second RL branch 
        vc1    #Voltage across capacitor in first RC branch 
        vc2    #Voltage across capacitor in second RC branch 
        η      #Integrator term for outer DC/DC PI controller
        κ      #Integrator term for inner DC/DC PI controller
        A      # First term for Pade approx
        B      # Second term for Pade approx
        C      # First term for Pade approx
        D      # Second term for Pade approx
        E      # Second term for Pade approx
        d_dc   # DC/DC converter duty cycle
        i_in   # Current flowing into the DC-link capacitor 
        v_conv # Magnitude of VSC output voltage

    end

    # Expressions
    pm = eg_d * ig_d + eg_q * ig_q  # AC Active Power Calculation
    qm = -eg_d * ig_q + eg_q * ig_d # AC Reactive Power Calculation
    ω_a = ωʳ + Dp * (pʳ - pf)  # Active Power Drop
    v_hat = vʳ + Dq * (qʳ - qf) # Reactive Power Drop
    v_iref_d = v_hat - rv * ig_d + ω_a * lv * ig_q # d-axis virtual impedance equation
    v_iref_q = -rv * ig_q - ω_a * lv * ig_d # q-axis virtual impedance equation
    i_hat_d = kvp * (v_iref_d - eg_d) + kvi * ξ_d - ω_a * cf * eg_q # Inner voltage controller d PI
    i_hat_q = kvp * (v_iref_q - eg_q) + kvi * ξ_q + ω_a * cf * eg_d # Inner voltage controller q PI
    v_md_bar = kip * (i_hat_d - is_d) + kii * γ_d - ω_a * lf * is_q # Inner current controller d PI
    v_mq_bar = kip * (i_hat_q - is_q) + kii * γ_q + ω_a * lf * is_d # Inner current controller q PI
    v_bar_mag = sqrt(v_md_bar^2 + v_mq_bar^2)
    v_md = (min(v_bar_mag, vdc) / v_bar_mag) * v_md_bar
    v_mq = (min(v_bar_mag, vdc) / v_bar_mag) * v_mq_bar
    p_inv = v_md * is_d + v_mq * is_q # Active power drawn from inverter
    q_inv = -v_md * is_q + v_mq * is_d # Reactive power drawn from inverter
    v_gd = (vl^2 / pl) * ig_d # d-qaxis Load voltage
    v_gq = (vl^2 / pl) * ig_q # q-axis load voltage
    i_ref = kpvb * (vdcʳ - vdc) + kivb * η # Current referecne for inner DC/DC PI controller 
    is_d_dot = (ωb / lf) * (v_md - eg_d) - (rf * ωb / lf) * is_d + ωb * ω_a * is_q # derivative of is_d
    is_q_dot = (ωb / lf) * (v_mq - eg_q) - (rf * ωb / lf) * is_q - ωb * ω_a * is_d # derivative of is_q
    i_pred = (Ts / vdc) * (v_md * is_d_dot + v_mq * is_q_dot) # Estiamtion of change in current flowing out of DC-link capacitor over DC/DC converter switching period

    model_rhs = [
        ### Grid forming equations
        #𝜕eg_d/𝜕t
        (ωb / cf) * (is_d - ig_d) + ω_a * ωb * eg_q
        #𝜕eg_q/𝜕t
        (ωb / cf) * (is_q - ig_q) - ω_a * ωb * eg_d
        #𝜕is_d/𝜕t
        (ωb / lf) * (v_md - eg_d) - (rf * ωb / lf) * is_d + ωb * ω_a * is_q
        #𝜕is_q/𝜕t
        (ωb / lf) * (v_mq - eg_q) - (rf * ωb / lf) * is_q - ωb * ω_a * is_d
        #𝜕ig_d/𝜕t
        (ωb / lt) * (eg_d - v_gd) - (rt * ωb / lt) * ig_d + ωb * ω_a * ig_q
        #𝜕ig_q/𝜕t
        (ωb / lt) * (eg_q - v_gq) - (rt * ωb / lt) * ig_q - ωb * ω_a * ig_d
        #𝜕pf/𝜕t
        ωz * (pm - pf)
        #𝜕qf/𝜕t
        ωz * (qm - qf)
        #𝜕ξ_d/𝜕t
        v_iref_d - eg_d
        #𝜕ξ_q/𝜕t
        v_iref_q - eg_q
        #𝜕γ_d/𝜕t
        i_hat_d - is_d
        #𝜕γ_q/𝜕t
        i_hat_q - is_q
        ### DC-side equations
        #∂vdc/∂t
        (ωb / cdc) * (i_in - p_inv / (vdc))
        #∂ibat/∂t
        (ωb / ldc) * (
            vb - vc1 - vc2 - rl2 * (ibat - il2) - rl1 * (ibat - il1) - rb0 * ibat -
            (1 - d_dc) * vdc
        )
        #∂il1/∂t
        (ωb / lb1) * (rl1 * (ibat - il1))
        #∂il2/∂t
        (ωb / lb2) * (rl2 * (ibat - il2))
        #∂vc1/∂t
        (ωb / cb1) * (ibat - vc1 / rc1)
        #∂vc2/∂t
        (ωb / cb2) * (ibat - vc2 / rc2)
        #∂η/∂t
        vdcʳ - vdc # Integrator for DC/DC outer PI controller
        #∂κ/dt
        i_ref + p_inv / (vdc) - i_in # Integrator for DC/DC inner PI controller
        # ∂A/dt
        (a1 / Ts) * A +
        (a2 / Ts^2) * B +
        kpib * (i_ref + p_inv / (vdc) - i_in) +
        kiib * κ +
        kpred * i_pred
        # ∂B/dt
        A # Second term in 2nd-order Pade approx.
        # ∂C/dt
        (b1 / Ts) * C +
        (b2 / Ts^2) * D +
        (b3 / Ts^3) * E +
        kpib * (i_ref + p_inv / (vdc) - i_in) +
        kiib * κ +
        kpred * i_pred
        # ∂D/dt
        C
        # ∂E/dt 
        D
        #Algebraic Eq.
        # DC/DC converter duty cycle
        -d_dc + min(
            0.9,
            (
                (1 / 2) * ((a2 / Ts) * A) +
                (1 / 2) * ((-2 * b1 / Ts) * C + (-2 * b3 / Ts^3) * E)
            ),
        )
        # Current flowing into DC-link capacitor
        -i_in + (vb * ibat - ibat^2 * rb0) / vdc
        # VSC AC terminal voltage 
        -v_conv + sqrt(v_md^2 + v_mq^2)
    ]

    # Temporary until SteadyState problems are resolved.
    model_lhs = [
        d(eg_d)
        d(eg_q)
        d(is_d)
        d(is_q)
        d(ig_d)
        d(ig_q)
        d(pf)
        d(qf)
        d(ξ_d)
        d(ξ_q)
        d(γ_d)
        d(γ_q)
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
        d(C)
        d(D)
        d(E)
        0
        0
        0
        #0
    ]

    return model_lhs, model_rhs, states, variables, params
end

function get_4th_order_dae_system()
    model_lhs, model_rhs, states, _, params = dae_model_4th_order(nothing)
    t = params[1]
    _eqs = model_lhs .~ model_rhs
    return MTK.ODESystem(_eqs, t, [states...], [params...][2:end])
end

function get_4th_order_nonlinear_system()
    _, model_rhs, _, variables, params = ode_model_4th_order(nothing)
    variable_count = length(variables)
    _eqs = zeros(length(model_rhs)) .~ model_rhs
    return MTK.NonlinearSystem(_eqs, [variables...], [params...][2:end])
end
