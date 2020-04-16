function ode_system!(dx,x,p,t)
    pm = x[1]*x[5] + x[2]*x[6]
    qm = -x[1]*x[6] + x[2]*x[5]
    ω_apc = p.ω_ref + p.Dp*(p.p_ref - x[7])
    v_hat = p.v_ref + p.Dq*(p.q_ref - x[7])
    v_iref_d = v_hat - p.rv*x[5] + ω_apc*p.lv*x[6]
    v_iref_q = - p.rv*x[6] - ω_apc*p.lv*x[5]
    i_hat_d = p.kvp*(v_iref_d -x[1]) + p.kvi*x[9] - ω_apc*p.cf*x[2]# + p.kif*x[5]
    i_hat_q = p.kvp*(v_iref_d -x[2]) + p.kvi*x[10] + ω_apc*p.cf*x[1]# + p.kif*x[6]
    v_md = p.kip*(i_hat_d -x[3]) + p.kii*x[11] - ω_apc*p.lf*x[4]# + p.kvf*x[1]
    v_mq = p.kip*(i_hat_q -x[4]) + p.kii*x[12] + ω_apc*p.lf*x[3]# + p.kvf*x[2]
    p_inv = v_md*x[3]+v_mq*x[4]
    q_inv = -v_md*x[4]+v_mq*x[3]
    v_gd = (p.vl^2/p.pl)*x[5]
    v_gq = (p.vl^2/p.pl)*x[6]
    ωg=ω_apc
    i_ref=p.kpvb*(p.vdc_ref-x[13])+p.kivb*x[15]
    i_in = (p.vb*x[14]-x[14]^2*p.req)/x[13]
    #d = p.kpib*(i_ref-i_in)+p.kiib*x[16]
    d = (-12/p.Ts)*x[17]+p.kpib*(i_ref-i_in)+p.kiib*x[16]
    #d = (-12/p.Ts)*x[17]+(72/p.Ts^2)*x[18]
    
### Grid forming equations
  #𝜕eg/𝜕t
    dx[1]= (p.ωb / p.cf)*(x[3] - x[5]) + ωg*p.ωb*x[2]
      #𝜕eq/𝜕t
    dx[2]= (p.ωb / p.cf)*(x[4] - x[6]) - ωg*p.ωb*x[1]
       #𝜕isd/𝜕t
    dx[3]= (p.ωb / p.lf)*(v_md - x[1]) - (p.rf*p.ωb/p.lf)*x[3] + p.ωb*ωg*x[4]
       #𝜕isq/𝜕t
    dx[4]= (p.ωb / p.lf)*(v_mq - x[2]) - (p.rf*p.ωb/p.lf)*x[4] - p.ωb*ωg*x[3]
         #𝜕igd/𝜕t
    dx[5]= (p.ωb / p.lt)*(x[1] - v_gd) - (p.rt*p.ωb/p.lt)*x[5] + p.ωb*ωg*x[6]
       #𝜕igq/𝜕t
    dx[6]= (p.ωb / p.lt)*(x[2] - v_gq) - (p.rt*p.ωb/p.lt)*x[6] - p.ωb*ωg*x[5] 
       #𝜕pf/𝜕t
    dx[7]=  p.ωz*(pm- x[7])
       #𝜕qf/𝜕t
    dx[8]= p.ωz*(qm - x[8])
      #𝜕δc/𝜕t
    #dx[9]= ωb*(ω_apc-ωg)
       #𝜕ξ_d/𝜕t
    dx[9]= v_iref_d - x[1]
       #𝜕ξ_q/𝜕t
    dx[10]= v_iref_q - x[2]
       #𝜕γ_d/𝜕t
    dx[11]=i_hat_d - x[3]
        #𝜕γ_q/𝜕t
    dx[12]= i_hat_q - x[4]
    
### DC-side equations
        #∂v_dc/∂t
    dx[13]=p.ωb*((i_in-p_inv/(2*x[13]))/(p.cdc));
        #∂i_batt/∂t
    dx[14]=(p.ωb/p.leq)*(p.vb-p.req*x[14]-(1-d)*x[13])    
        #∂η/∂t Integrator for DC/DC outer PI controller
    dx[15]=p.vdc_ref-x[13]
        #∂κ/dt Integrator for DC/DC inner PI controller
    dx[16]=i_ref-i_in
        # First term in Pade approximation
    dx[17]=(-6/p.Ts)*x[17] + (-12/p.Ts^2)*x[18] + p.kpib*(i_ref-i_in)+p.kiib*x[16]
        # Second term in Pade approx.
    dx[18]=x[17]
end
