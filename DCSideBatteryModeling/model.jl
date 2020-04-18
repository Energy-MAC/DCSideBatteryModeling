MTK.@parameters begin
    t
    ω_ref
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
pm = eg_d*ig_d + eg_q*ig_q
qm = -eg_d*ig_q + eg_q*ig_d
ω_a = ω_ref + p.Dp*(p.p_ref - x[7])
v_hat = p.v_ref + p.Dq*(p.q_ref - x[7])
v_iref_d = v_hat - p.rv*ig_d(t) + ω_a*p.lv*x[6]
v_iref_q = - p.rv*x[6] - ω_a*p.lv*ig_d(t)
i_hat_d = p.kvp*(v_iref_d -eg_d(t)) + p.kvi*x[9] - ω_a*p.cf*eg_q(t)# + kif*ig_d(t)
i_hat_q = p.kvp*(v_iref_d -eg_q(t)) + p.kvi*x[10] + ω_a*p.cf*eg_d(t)# + kif*x[6]
v_md = p.kip*(i_hat_d -x[3]) + p.kii*x[11] - ω_a*p.lf*x[4]# + kvf*eg_d(t)
v_mq = p.kip*(i_hat_q -x[4]) + p.kii*x[12] + ω_a*p.lf*x[3]# + kvf*eg_q(t)
inv = v_md*x[3]+v_mq*x[4]
q_inv = -v_md*x[4]+v_mq*x[3]
v_gd = (p.vl^2/p.pl)*ig_d(t)
v_gq = (p.vl^2/p.pl)*x[6]
ωg=ω_a
i_ref=p.kpvb*(p.vdc_ref-x[13])+p.kivb*x[15]
i_in = (p.vb*x[14]-x[14]^2*p.req)/x[13]
d = (-12/p.Ts)*x[17]+p.kpib*(i_ref-i_in)+p.kiib*x[16]

model = [



]
