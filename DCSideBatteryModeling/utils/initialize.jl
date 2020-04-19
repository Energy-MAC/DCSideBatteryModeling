
function get_ini_sys()
    model, states, params = _get_system()
    ini_sys = MTK.NonlinearSystem(model, [states...], [params...])
end

#=
function get_IC()
    #Initial Conditions

    1. Solve PowerFlow with PSY

    initial = [
        1.0, #egd1
        0.0, #egq1
        0.5, #isd1
        0, #isq1
        0.5, #igd1
        0,#igq1
        0.5,
        0.0,
        0.0, #ξ_d1
        0.0, #ξ_q1
        0.0, #γ_d1
        0.0, #γ_q1
        1.0, #v_dc
        0.5, #i_batt
        0.0, #η
        0.0, #κ
        0.0,
        0.0,
    ]

    #= From MTK documentation
nlsys_func = generate_function(ns)[2] # second is the inplace version

f = eval(nlsys_func)
du = zeros(3); u = ones(3)
params = (10.0,26.0,2.33)
f(du,u,params)
du

j_func = generate_jacobian(ns)[2] # second is in-place
j! = eval(j_func)

using NLsolve
nlsolve((out, x) -> f(out, x, params), (out, x) -> j!(out, x, params), ones(3))
    #=
       sys_solve = NLsolve.nlsolve(system!, initial, method = :trust_region) #Solve using initial guess x0
       #print(sys_solve)
       return sys_solve
    end
=#
