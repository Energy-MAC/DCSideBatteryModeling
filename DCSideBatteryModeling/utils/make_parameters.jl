
Ub = 690;
fb = 60;
ωb = 2 * pi * fb;
Sb = 0.8e6;
Vb = 690;

function get_params(Ub, fb, ωb, Sb, Vb)
    # System base for AC side
    Ib = Sb / Ub
    Zb = Vb / Ib
    Lb = Zb / ωb
    Cb = 1 / (Zb * ωb)

    # System base for DC side
    Vb_dc = 2 * Vb
    Ib_dc = Sb / Vb_dc
    Zb_dc = Vb_dc / Ib_dc
    Lb_dc = Zb_dc / ωb
    Cb_dc = 1 / (Zb_dc * ωb)

    # Grid impadance (Currently not used)
    lg = 0.2
    rg = 0.01

    #Reference set-point inputs
    p_ref = 0.5
    q_ref = 0
    v_ref = 1
    ω_ref = 1
    vdc_ref = 2.02 * 690 / Vb_dc

    # Load at rated voltage
    vl = 1
    pl = 0.5

    # Filter parameters
    lf = 0.08
    cf = 0.074
    rf = 0.003

    # Filterwing frequency
    ωz = 0.1 * ωb

    # Transformer Paramters
    rt = 0.01
    lt = 0.2

    # Outer Control Loops
    # Active Power Droop Control
    Dp = 0.02

    # Reactive Power Droop Control
    Dq = 0.001

    # Inner Control Loops
    # SRF Current Control
    kip = 1.27
    kii = 14.3
    kffi = 0.0

    # SRF Voltage Control
    kvp = 0.59
    kvi = 736
    kffv = 1.0

    # Virtual Impedance
    rv = 0
    lv = 0.2

    # DC source Parameters
    leq = 5e-3 / Lb_dc
    req = 0.0 / Zb_dc
    vb = 1000 / Vb_dc
    cdc = 100e-6 / Cb_dc

    # DC/DC converter controller parameters
    kpvb = 0.6
    kivb = 4
    kpib = 0.3863
    kiib = 10.34
    Ts = 1 / (3.2e3)

    return get_params(
        ωb,
        lg,
        rg,
        p_ref,
        q_ref,
        v_ref,
        ω_ref,
        vl,
        pl,
        lf,
        cf,
        rf,
        ωz,
        rt,
        lt,
        Dp,
        Dq,
        kip,
        kii,
        kffi,
        kvp,
        kvi,
        kffv,
        rv,
        lv,
        leq,
        req,
        vb,
        vdc_ref,
        cdc,
        kivb,
        kpvb,
        kpib,
        kiib,
        Ts,
    )
end
