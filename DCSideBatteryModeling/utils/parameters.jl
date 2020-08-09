function instantiate_parameters(
    system::PSY.System,
    model = get_4th_order_nonlinear_system(),
)
    # Base quantities
    ωb,      # Base Frequency
    Cb_dc,   # DC-side base capacitance
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
    kip,     # VSC Current control propotional gain
    kii,     # VSC Current control integral gain
    kffi,    # VSC Current control feed-forward gain
    # SRF Voltage Control
    kvp,     # VSC Voltage control propotional gain
    kvi,     # VSC Voltage control integral gain
    kffv,    # VSC Voltage control feed-forward  gain
    # Virtual Impedance
    rv,      # Virutal resistance
    lv,      # Virutal inductance
    # DC Source Parameters
    ldc,     # DC/DC inductance
    req,     # Equivalent resistance for 0th-order model
    rb0,     # Battery steady-state resistance
    lb1,     # Inductacne of 1st RL branch in battery model
    rl1,     # Resistance of 1st RL branch in battery model
    lb2,     # Inductacne of 2nd RL branch in battery model
    rl2,     # Resistance of 2nd RL branch in battery model
    cb1,     # Capacitance of first RC branch
    rc1,     # Resistance of first RC branch
    cb2,     # Capacitance of second RC branch
    rc2,     # Resistance of second RC branch
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

    Ub = 480 # AC phase base voltage
    fb = 60 # Base frequency
    _ωb = 2 * pi * fb # Base angular frequency
    Sb = 0.2e6 # Base power
    Vb = 480 # AC phase base voltage
    # System base for AC side
    Ib = Sb / Ub # AC Base current
    Zb = Vb / Ib # AC Base impedance
    Lb = Zb / _ωb # AC Base inductance
    Cb = 1 / (Zb * _ωb) # AC Base capacitance
    # System base for DC side
    Vb_dc = (sqrt(8) / sqrt(3)) * Vb # DC-side base voltage
    Ib_dc = Sb / Vb_dc # DC-side base current
    Zb_dc = Vb_dc / Ib_dc # DC-side base impedance
    Lb_dc = Zb_dc / _ωb # DC-side base inductance
    _Cb_dc = 1 / (Zb_dc * _ωb) # DC-side base capacitance

    p = [
        ωb => _ωb # System base angular frequency
        Cb_dc => _Cb_dc # DC-side base capacitance
        # Grid impadance (Currently not used)
        lg => 0.2 # Grid inductance (p.u.)
        rg => 0.01 # Grid resistance (p.u.)
        #Reference set-point inputs
        pʳ => 0.2 # VSC active power reference set-point (p.u.)
        qʳ => 0 # VSC reactive power reference set-point (p.u.)
        vʳ => 1.01 # VSC voltage reference set-point (p.u.)
        ωʳ => 1 # VSC Reference frequency (p.u.)
        # Load
        vl => 1 # Rated voltage (p.u.)
        pl => 0.2 # Active power at rated voltage (p.u.)
        # Filter parameters
        lf => 0.08  # VSC Filter inductance (p.u.)
        cf => 0.074 # VSC Filter capacitance (p.u.)
        rf => 0.003 # VSC filer resistance (p.u.)
        # Filtering frequency
        ωz => 0.1 * _ωb # (p.u.)
        # Transformer Parameters
        rt => 0.01 # Transformer resistance
        lt => 0.2  # Transfer inductance
        # Outer Control Loops
        Dp => 0.02 # VSC active power droop gain
        Dq => 0.001 # VSC reactive power droop gain
        # SRF Current Control
        kip => 1.27 # VSC current control propotional gain
        kii => 14.3 # VSC current control integral gain
        kffi => 0.0 # VSC current control feed-forward gain
        # SRF Voltage Control
        kvp => 0.59 # VSC Voltage control propotional gain
        kvi => 736  # VSC Voltage control integral gain
        kffv => 1.0 # VSC Voltage control feed-forward gain
        # Virtual Impedance
        rv => 0     # Virtual resistance in p.u.
        lv => 0.2   # Virtual inductance in p.u.
        # DC source Parameters
        ldc => 3e-3 / Lb_dc # DC/DC inductor in p.u.
        req => (1.5e-3 + 2.2e-3 + 0.55e-3) / Zb_dc # Equivalent battery resistance for 0th-order model in p.u.
        rb0 => 1.5e-3 / Zb_dc  # Battery steady-state resistance
        lb1 => 35e-9 / Lb_dc # Inductacne of 1st RL branch in battery model in p.u.
        rl1 => 95e-3 / Zb_dc # Resistance of 1st RL branch in battery model in p.u.
        lb2 => 15e-9 / Lb_dc # Inductacne of 2nd RL branch in battery model in p.u.
        rl2 => 0.4e-3 / Zb_dc # Resistance of 2nd RL branch in battery model in p.u.
        cb1 => 0.55 / _Cb_dc # Capacitance of first RC branch in p.u.
        rc1 => 2.2e-3 / Zb_dc # Resistance of first RC branch in p.u.
        cb2 => 22700 / _Cb_dc # Capacitance of second RC branch in p.u. 
        rc2 => 0.55e-3 / Zb_dc # Resistance of second RC branch in p.u.
        vb => 370 / Vb_dc  # Battery voltage in p.u.
        cdc => 2000e-6 / _Cb_dc # DC-link capacitance in p.u.
        # DC/DC converter controller parameters
        vdcʳ => 850 / Vb_dc # DC voltage reference
        kpvb => 0.6 # Outer-loop DC/DC proportional gain
        kivb => 4   # Outer-loop DC/DC integral gain
        kpib => 0.3863 # Inner-loop DC/DC Controller proportional gain
        kiib => 10.34 # Inner-loop DC/DC Controller integral gain
        kpred => 0 #Gain on one step predictor
        a1 => -6    # First coefficient of 2nd order Pade approximation
        a2 => -12    # Second co-efficient of 2nd order Pade approxmiation
        b1 => -12    # First coefficient of 3rd order Pade approximation
        b2 => -60    # Second co-efficient of 3rd order Pade approxmiation
        b3 => -120    # Third coefficient of 3rd order Pade approximation
        Ts => 1 / (3.2e3) # Get using PSY
    ]
    return p
end
