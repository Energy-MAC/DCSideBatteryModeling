function get_system()

    # Model Parameters
    MTK.@parameters begin
        t
        ωg
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
        leq
        req
        vb
        cdc
        # DC/DC converter controller parameters
        vdcʳ    # DC Voltage reference
        kivb
        kpvb
        kpib
        kiib
        Ts # Sampling time
    end

    MTK.@derivatives d'~t

    # Definition of the states
    MTK.@variables begin
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
        η(t)    #Integrator term for outer DC/DC PI controller
        κ(t)    #Integrator term for inner DC/DC PI controller
        # TODO: Verify in the nomenclature equation is the appropiate for each term of the Pade approximation
        M(t)    # First term for Pade approx
        L(t)    # Second term for Pade approx
    end

    # Expressions
    pm = eg_d * ig_d + eg_q * ig_q  # AC Active Power Calculation
    qm = -eg_d * ig_q + eg_q * ig_d # AC Reactive Power Calculation
    ω_a = ωʳ + Dp * (pʳ - pf)  # Active Power Drop
    # TODO: Original model had pf here. Verify
    v_hat = vʳ + Dq * (qʳ - qf) # Reactive Power Drop
    v_iref_d = v_hat - rv * ig_d + ω_a * lv * pf # Inner voltage controller d PI
    v_iref_q = -rv * pf - ω_a * lv * ig_d # Inner voltage controller q PI
    i_hat_d = kvp * (v_iref_d - eg_d) + kvi * ξ_d - ω_a * cf * eg_q # Inner current controller d PI
    i_hat_q = kvp * (v_iref_d - eg_q) + kvi * ξ_q + ω_a * cf * eg_d # Inner current controller q PI
    v_md = kip * (i_hat_d - is_d) + kii * γ_d - ω_a * lf * is_q
    v_mq = kip * (i_hat_q - is_q) + kii * γ_q + ω_a * lf * is_d
    p_inv = v_md * is_d + v_mq *is_q
    q_inv = -v_md * is_q + v_mq * is_d
    v_gd = (vl^2 / pl) * ig_d
    v_gq = (vl^2 / pl) * ig_q
    ωg = ω_a
    i_ref = kpvb * (vdcʳ - vdc) + kivb * η
    i_in = (vb * ibat - ibat^2 * req) / vdc
    d_dc = (-12 / Ts) * M + kpib * (i_ref - i_in) + kiib * κ

    model = [
        ### Grid forming equations
        #𝜕eg_d/𝜕t
        d(eg_d) ~ (ωb / cf) * (is_d - ig_d) + ωg * ωb * eg_q
        #𝜕eg_q/𝜕t
        d(eg_q) ~ (ωb / cf) * (is_q - ig_q) - ωg * ωb * eg_d
        #𝜕is_d/𝜕t
        d(is_d) ~ (ωb / lf) * (v_md - eg_d) - (rf * ωb / lf) * is_d + ωb * ωg * is_q
        #𝜕is_q/𝜕t
        d(is_q) ~ (ωb / lf) * (v_mq - eg_q) - (rf * ωb / lf) * is_q - ωb * ωg * is_d
        #𝜕ig_d/𝜕t
        d(ig_d) ~ (ωb / lt) * (eg_d - v_gd) - (rt * ωb / lt) * ig_d + ωb * ωg * ig_q
        #𝜕ig_q/𝜕t
        d(ig_q) ~ (ωb / lt) * (eg_q - v_gq) - (rt * ωb / lt) * ig_q - ωb * ωg * ig_d
        #𝜕pf/𝜕t
        d(pf) ~ ωz * (pm - pf)
        #𝜕qf/𝜕t
        d(qf) ~ ωz * (qm - qf)
        #𝜕ξ_d/𝜕t
        d(ξ_d) ~ v_iref_d - eg_d
        #𝜕ξ_q/𝜕t
        d(ξ_q) ~ v_iref_q - eg_q
        #𝜕γ_d/𝜕t
        d(γ_d) ~ i_hat_d - is_d
        #𝜕γ_q/𝜕t
        d(γ_q) ~ i_hat_q - is_q
        ### DC-side equations
        #∂vdc/∂t
        d(vdc) ~ ωb * ((i_in - p_inv / (2 * vdc)) / (cdc))
        #∂ibat/∂t
        d(ibat) ~ (ωb / leq) * (vb - req * ibat - (1 - d_dc) * vdc)
        #∂η/∂t
        d(η) ~ vdcʳ - vdc # Integrator for DC/DC outer PI controller
        #∂κ/dt
        d(κ) ~ i_ref - i_in # Integrator for DC/DC inner PI controller
        # ∂M/dt
        d(M) ~
            (-6 / Ts) * M +
            (-12 / Ts^2) * L +
            kpib * (i_ref - i_in) +
            kiib * ibat # First term in Pade approximation
        # ∂M/dt
        d(L) ~ M # Second term in Pade approx.
    ]

    return MTK.ODESystem(model)
end
