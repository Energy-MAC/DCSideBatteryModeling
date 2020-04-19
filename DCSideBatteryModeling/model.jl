MTK.@parameters begin
    t
    # AC side quantities
    ωb  # Base Frequency
    lg  # Grid reactance
    rg  # Grid resistance
    pʳ  # Active Power Setpoint
    qʳ  # Reactive Power Setpoint
    vʳ  # Voltage Setpoint
    ωʳ  # Frequency Setpoint
    vl
    pl
    lf
    cf
    rf
    ωz
    rt
    lt
    Dp  # Active Power Damping
    Dq  # Reactive Power Damping
    kip
    kii
    kffi
    kvp
    kvi
    kffv
    rv
    lv
    # DC Side quantities
    leq
    req
    vb
    vdcʳ
    cdc
    kivb
    kpvb
    kpib
    kiib
    Ts # Sampling time
end

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

MTK.@derivatives d'~t

# Expressions
pm = eg_d * ig_d + eg_q * ig_q  # AC Active Power Calculation
qm = -eg_d * ig_q + eg_q * ig_d # AC Reactive Power Calculation
ω_a = ωʳ + Dp * (pʳ - pf)  # Active Power Drop
# TODO: Original model had x[7] here. Verify
v_hat = vʳ + Dq * (qʳ - qf) # Reactive Power Drop
v_iref_d = v_hat - rv * ig_d + ω_a * lv * pf # Inner voltage controller d PI
v_iref_q = -rv * pf - ω_a * lv * ig_d # Inner voltage controller q PI
i_hat_d = kvp * (v_iref_d - eg_d) + kvi * ξ_d - ω_a * cf * eg_q # Inner current controller d PI
i_hat_q = kvp * (v_iref_d - eg_q) + kvi * ξ_q + ω_a * cf * eg_d # Inner current controller q PI
v_md = kip * (i_hat_d - x[3]) + kii * x[11] - ω_a * lf * x[4]
v_mq = kip * (i_hat_q - x[4]) + kii * x[12] + ω_a * lf * x[3]
inv = v_md * x[3] + v_mq * x[4]
q_inv = -v_md * x[4] + v_mq * x[3]
v_gd = (vl^2 / pl) * ig_d(t)
v_gq = (vl^2 / pl) * x[6]
ωg = ω_a
i_ref = kpvb * (vdc_ref - x[13]) + kivb * x[15]
i_in = (vb * x[14] - x[14]^2 * req) / x[13]
d = (-12 / Ts) * x[17] + kpib * (i_ref - i_in) + kiib * x[16]

model = []
