function instantiate_parameters(system::PSY.System, model = get_4th_order_nonlinear_system())
    # TODO: Automate better with PSY getter functions
    # AC side quantities
    ωb,      # Base Frequency
    # Grid impadance
    lg,      # Grid reactance
    rg,      # Grid resistance
    #Reference set-point input
    pʳ,      # Active Power Setpoint
    qʳ,      # Reactive Power Setpoint
    vʳ,      # Voltage Setpoint
    ωʳ,      # Frequency Setpoint
    # Load at rated voltage
    vl,      # Load rated voltage
    pl,      # Load rated power
    # Filter parameters
    lf,      # Filter reactance
    cf,      # Filter capacitance
    rf,      # Filter Resistance
    # Filtering frequency
    ωz,      # Filtering frequency
    # Transformer Parameters
    rt,      # Transformer resistance
    lt,      # Transformer reactance
    # OuterControl Loops
    Dp,      # Active Power Droop
    Dq,      # Reactive Power Droop
    # SRF Current Control
    kip,     # Current control propotional gain
    kii,     # Current control integral gain
    kffi,    # Current control feed-forward gain
    # SRF Voltage Control
    kvp,     # Voltage control propotional gain
    kvi,     # Voltage control integral gain
    kffv,    # Voltage control feed-forward  gain
    # Virtual Impedance
    rv,
    lv,
    # DC Source Parameters
    ldc,     # DC/DC inductance
    req,     #
    rb0,     # Battery steady-state resistance
    lb1,     # Inductacne of 1st RL branch in battery model
    rl1,     # Resistance of 1st RL branch in battery model
    lb2,     # Inductacne of 2nd RL branch in battery model
    rl2,     # Resistance of 2nd RL branch in battery model
    cb1,     #
    rc1,     #
    cb2,     #
    rc2,     #
    vb,      # Battery Voltage
    cdc,     # Dc-side capacitance
    # DC/DC converter controller parameters
    vdcʳ,    # DC Voltage reference
    kpvb,    # DC/DC Voltage control integral gain
    kivb,    # DC/DC Voltage control propotional gain
    kpib,    # DC/DC Current control propotional gain
    kiib,    # DC/DC Current control Integral gain
    kpred,   # Gain on one-step predictor
    a1,      # First coefficient of 2nd order Pade approximation
    a2,      # Second co-efficient of 2nd order Pade approxmiation
    b1,      # First coefficient of 3rd order Pade approximation
    b2,      # Second coefficient of 3rd order Pade approximation
    b3,      # Third coefficient of 3rd order Pade approximation   
    Ts = MTK.parameters(model) ## DC/DC controller time delay

    Ub = 480 # Get using PSY
    fb = 60 # Get using PSY.
    _ωb = 2 * pi * fb
    Sb = 0.2e6 # Get using PSY
    Vb = 480 # Get using PSY
    # System base for AC side
    Ib = Sb / Ub
    Zb = Vb / Ib
    Lb = Zb / _ωb
    Cb = 1 / (Zb * _ωb)
    # System base for DC side
    Vb_dc = (sqrt(8)/sqrt(3))*Vb
    Ib_dc = Sb / Vb_dc
    Zb_dc = Vb_dc / Ib_dc
    Lb_dc = Zb_dc / _ωb
    Cb_dc = 1 / (Zb_dc * _ωb)

    p = [
        ωb => _ωb
        # Grid impadance (Currently not used)
        lg => 0.2 # Get using PSY
        rg => 0.01 # Get using PSY
        #Reference set-point inputs
        pʳ => 0.5 # Get using PSY
        qʳ => 0 # Get using PSY
        vʳ => 1.01 # Get using PSY
        ωʳ => 1 # Reference frequency
        # Load at rated voltage
        vl => 1 # Get using PSY
        pl => 0.2 # Get using PSY
        # Filter parameters
        lf => 0.08  # Get using PSY
        cf => 0.074 # Get using PSY
        rf => 0.003 # Get using PSY
        # Filtering frequency
        ωz => 0.1 * _ωb
        # Transformer Parameters
        rt => 0.01 # Get using PSY
        lt => 0.2  # Get using PSY
        # Outer Control Loops
        Dp => 0.02 # Get using PSY
        Dq => 0.001 # Get using PSY
        # SRF Current Control
        kip => 1.27 # Get using PSY
        kii => 14.3 # Get using PSY
        kffi => 0.0 # Get using PSY
        # SRF Voltage Control
        kvp => 0.59 # Get using PSY
        kvi => 736  # Get using PSY
        kffv => 1.0 # Get using PSY
        # Virtual Impedance
        rv => 0     # Get using PSY
        lv => 0.2   # Get using PSY
        # DC source Parameters
        ldc => 3e-3 / Lb_dc # Get using PSY
        req => (1.5e-3+2.2e-3+0.55e-3) / Zb_dc
        rb0 => 1.5e-3 / Zb_dc  # Battery steady-state resistance
        lb1 => 35e-9 / Lb_dc # Inductacne of 1st RL branch in battery model
        rl1 => 95e-3 / Zb_dc # Resistance of 1st RL branch in battery model
        lb2 => 15e-9 / Lb_dc # Inductacne of 2nd RL branch in battery model
        rl2 => 0.4e-3 / Zb_dc # Resistance of 2nd RL branch in battery model
        cb1 => 0.55/Cb_dc #
        rc1 => 2.2e-3/Zb_dc #
        cb2 => 22700/Cb_dc #
        rc2 => 0.55e-3/Zb_dc #
        vb => 370 / Vb_dc  # Get using PSY
        cdc => 4000e-6 / Cb_dc # Get using PSY
        # DC/DC converter controller parameters
        vdcʳ => 850/ Vb_dc # Get using PSY
        kpvb => 0.6 # Get using PSY
        kivb => 4   # Get using PSY
        kpib => 0.3863 # Get using PSY
        kiib => 10.34 # Get using PSY
        kpred => 0.01 #Gain on one step predictor
        a1 => -6    # First coefficient of 2nd order Pade approximation
        a2 => -12    # Second co-efficient of 2nd order Pade approxmiation
        b1 => -12    # First coefficient of 3rd order Pade approximation
        b2 => -60    # Second co-efficient of 3rd order Pade approxmiation
        b3 => -120    # Third coefficient of 3rd order Pade approximation
        Ts => 1 / (3.2e3) # Get using PSY
    ]
    return p
end
