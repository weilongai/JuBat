function sP2D(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    variables = sP2D_variables(case, yt, t)
    param = case.param
    mesh_np = case.mesh["negative particle"]
    mesh_pp = case.mesh["positive particle"]
    mesh_el = case.mesh["electrolyte"]
    if jacobi == "constant" && param.NE.M_d != [] # no need to update M and K
        M_ne_d = param.NE.M_d
        K_ne_d = param.NE.K_d
        M_pe_d = param.PE.M_d
        K_pe_d = param.PE.K_d
    else
        M_ne_d, K_ne_d = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen)
        M_pe_d, K_pe_d = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen)
        M_ne_d = M_ne_d .* param.scale.ts_n / param_dim.scale.t0
        M_pe_d = M_pe_d .* param.scale.ts_p / param_dim.scale.t0   
    end
    M_el_d, K_el_d = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, variables)   
    M_el_d = M_el_d .* param.scale.te / param_dim.scale.t0  

    # # need to update source term
    phi_new, variables = sP2D_potentials(case, yt, t, variables)
    F = sP2D_mass_BC(case, variables)
    M = blockdiag(M_ne_d, M_pe_d, M_el_d)
    K = blockdiag(K_ne_d, K_pe_d, K_el_d)
    return M, K, F, variables, phi_new
end

function sP2D_mass_BC(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param

    # electrode particle diffusion source term
    j_n = variables["negative electrode interfacial current density"]
    j_p = variables["positive electrode interfacial current density"]
    flux_np = zeros(Float64, case.mesh["negative particle"].nlen, 1)
    v_np_surf = case.index["negative particle surface lithium concentration"]
    flux_np[v_np_surf] = - j_n * param.NE.rs^2
    flux_pp = zeros(Float64, case.mesh["positive particle"].nlen, 1)
    v_pp_surf = case.index["positive particle surface lithium concentration"] .- case.mesh["negative particle"].nlen
    flux_pp[v_pp_surf] = - j_p * param.PE.rs^2

    # electrolyte diffusion source term
    mesh_el = case.mesh["electrolyte"]
    n_sp_gs = case.opt.Ns * mesh_el.gs.order
    j_n_gs = variables["negative electrode interfacial current at Gauss point"]
    j_p_gs = variables["positive electrode interfacial current at Gauss point"]
    aj_el_gs = [j_n_gs * param.NE.as;  zeros(n_sp_gs, 1); j_p_gs * param.PE.as] 
    coeff_el = mesh_el.gs.weight .* mesh_el.gs.detJ .* aj_el_gs .* (1 - param.EL.tplus)
    Vi_el = mesh_el.element[mesh_el.gs.ele,:]
    flux_el = Assemble1D(Vi_el, mesh_el.gs.Ni, coeff_el, mesh_el.nlen)

    # assemble all fluxes
    return [flux_np; flux_pp; flux_el]
end

function sP2D_potentials(case::Case, yt::Array{Float64}, t::Float64, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param  
    mesh_np = case.mesh["negative particle"]
    mesh_pp = case.mesh["positive particle"]
    mesh_el = case.mesh["electrolyte"]
    mesh_sp = case.mesh["separator"]
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    I_app =case.opt.Current(t * param.scale.t0) / param.scale.I_typ
    j0_n_gs = variables["negative electrode exchange current density at Gauss point"]
    j0_p_gs = variables["positive electrode exchange current density at Gauss point"]
    u_n_gs = variables["negative electrode open circuit potential at Gauss point"]
    u_p_gs = variables["positive electrode open circuit potential at Gauss point"]
    T = variables["temperature"]  
    ce_n = variables["electrolyte lithium concentration in negative electrode"]
    ce_p = variables["electrolyte lithium concentration in positive electrode"]
    ce_sp = variables["electrolyte lithium concentration in separator"] 

    phis_n = 0
    phis_p0 = u_p_gs[end]
    xn_gs = mesh_ne.gs.x 
    xp_gs = mesh_pe.gs.x 
    Li = [param.NE.thickness, param.SP.thickness, param.PE.thickness]
    dcedx = (ce_n[end] - ce_n[end - 1]) / (mesh_ne.node[end] - mesh_ne.node[end - 1])
    ce_v = (ce_n[end] + ce_n[end-1] ) / 2
    kn = 2 * T * (1 - param.EL.tplus) * param.EL.dlnf_dlnc(ce_v) / ce_v * dcedx
    kappa_n_eff = param.EL.kappa(ce_v) * param.NE.eps ^ param.NE.brugg
    kn -= I_app / kappa_n_eff
    dcedx = (ce_p[2] - ce_p[1]) / (mesh_pe.node[2] - mesh_pe.node[1])
    ce_v = ( ce_p[1] + ce_p[2]) / 2
    kp = 2 * T * (1 - param.EL.tplus) * param.EL.dlnf_dlnc(ce_v) / ce_v * dcedx
    kappa_p_eff = param.EL.kappa(ce_v) * param.PE.eps ^ param.PE.brugg
    kp -= I_app / kappa_p_eff
    v_sp = Int64((mesh_sp.nlen + 1) / 2)
    dcedx = (ce_sp[v_sp + 1] - ce_sp[v_sp - 1]) / (mesh_sp.node[v_sp + 1] - mesh_sp.node[v_sp - 1])
    ks = (kn * param.NE.eps ^ param.NE.brugg + kp * param.PE.eps ^ param.PE.brugg) / 2 / param.SP.eps ^ param.SP.brugg
    ki = [kn, ks, kp]
    phie_n_gs_rel = phie_fit(xn_gs, Li, ki, 0.0)
    phie_p_gs_rel = phie_fit(xp_gs, Li, ki, 0.0)
    eta_n_gs_rel = phis_n .- phie_n_gs_rel - u_n_gs
    eta_p_gs_rel = phis_p0 .- phie_p_gs_rel - u_p_gs
    I_np = IntV(param.NE.as .* j0_n_gs .* exp.(0.5 * eta_n_gs_rel / T), mesh_ne)
    I_nn = IntV(param.NE.as .* j0_n_gs .* exp.(-0.5 * eta_n_gs_rel / T), mesh_ne)
    I_pp = IntV(param.PE.as .* j0_p_gs .* exp.(0.5 * eta_p_gs_rel / T), mesh_pe)
    I_pn = IntV(param.PE.as .* j0_p_gs .* exp.(-0.5 * eta_p_gs_rel / T), mesh_pe)
    phie0 = - 2.0 * T * log((I_app + sqrt(4.0 * I_np * I_nn + I_app^2.0)) / 2.0 / I_np)
    phis_p = 2.0 * T * log((- I_app + sqrt(4.0 * I_pp * I_pn + I_app^2.0)) / 2.0 / I_pp) + phie0 + phis_p0
    phie = phie_fit(mesh_el.node, Li, ki, phie0)
    phi_new = [phis_n * ones(mesh_ne.nlen,1); phis_p * ones(mesh_pe.nlen,1); phie]
    yt_new = [yt[1:mesh_np.nlen + mesh_pp.nlen + mesh_el.nlen,1]; phi_new]
    variables = sP2D_variables(case, yt_new, t)
    return phi_new, variables
end

function phie_fit(x::Union{Array{Float64},Float64}, Li::Array{Float64}, ki::Array{Float64}, phie0::Float64, typef::String="sin")
    # a function to fit potential curve in electrolyte 
    Ln, Ls, Lp = Li
    kn, ks, kp = ki
    L = Ln + Ls + Lp
    phie = zeros(length(x))
    # # # one-stage fitting - quad function
    if typef == "sin"
        for i in eachindex(x)
            if x[i] <= Ln
                phie[i]= 0.5 * kn / Ln * x[i]^2 - kn * Ln/2 - ks * Ls/2 + phie0
            elseif x[i] < Ln + Ls
                phie[i]= ks * (x[i] - Ln - Ls/2) + phie0
            else
                phie[i]= -0.5 * kp / Lp * (x[i] - L)^2 + kp * Lp/2 + ks * Ls/2 + phie0
            end
        end
    elseif  typef == "quad"
    # one-stage fitting - sin function
        pi = 3.1415
        for i in eachindex(x)
            if x[i] <= Ln
                phie[i]= - 2 * kn * Ln / pi * cos(x[i] * pi / 2 / Ln) - ks * Ls/2 + phie0
            elseif x[i] < Ln + Ls
                phie[i]= ks * (x[i] - Ln - Ls/2) + phie0
            else
                phie[i]= 2 * kp * Lp /pi * sin((x[i] - Ln - Ls) / Lp * pi / 2) + ks * Ls/2 + phie0
            end
        end
    end

    return phie
end

function sP2D_variables(case::Case, yt::Array{Float64}, t::Float64)
    param = case.param
    variables = StandardVariables(case, 1)
    var_list = collect(keys(case.index))
    for i in var_list
        variables[i] = yt[case.index[i]] # hcat converts vector to matrix
    end
    if "temperature" in var_list
        T = yt[case.index["temperature"]]
    else
        T = case.param.cell.T0
    end
    gs_ne = case.mesh["negative electrode"].gs
    gs_pe = case.mesh["positive electrode"].gs
    gs_sp = case.mesh["separator"].gs
    element_ne = case.mesh["negative electrode"].element
    element_pe = case.mesh["positive electrode"].element
    element_sp = case.mesh["separator"].element

    csn_surf = variables["negative particle surface lithium concentration"]
    csp_surf = variables["positive particle surface lithium concentration"]
    phis_n = variables["negative electrode potential"]
    phis_p = variables["positive electrode potential"]
    phie_n = variables["electrolyte potential in negative electrode"]
    phie_p = variables["electrolyte potential in positive electrode"]
    ce_n = variables["electrolyte lithium concentration in negative electrode"]
    ce_p = variables["electrolyte lithium concentration in positive electrode"]
    ce_sp = variables["electrolyte lithium concentration in separator"]
    u_n = param.NE.U(csn_surf)
    u_p = param.PE.U(csp_surf)
    eta_p = phis_p - phie_p - u_p
    eta_n = phis_n - phie_n - u_n
    j0_n =  param.NE.k * Arrhenius(param.NE.Eac_k, T) .* abs.(csn_surf .* (1.0 .- csn_surf) .* ce_n) .^ 0.5
    j0_p =  param.PE.k * Arrhenius(param.PE.Eac_k, T) .* abs.(csp_surf .* (1.0 .- csp_surf) .* ce_p) .^ 0.5
    j_n = j0_n .* sinh.(0.5 .* eta_n / T) * 2.0
    j_p = j0_p .* sinh.(0.5 .* eta_p / T) * 2.0

    csn_surf_gs = sum(gs_ne.Ni .* csn_surf[element_ne[gs_ne.ele,:]], dims=2)
    csp_surf_gs = sum(gs_pe.Ni .* csp_surf[element_pe[gs_pe.ele,:]], dims=2)
    phie_n_gs = sum(gs_ne.Ni .* phie_n[element_ne[gs_ne.ele,:]], dims=2)
    phie_p_gs = sum(gs_pe.Ni .* phie_p[element_pe[gs_pe.ele,:]], dims=2)
    ce_n_gs = sum(gs_ne.Ni .* ce_n[element_ne[gs_ne.ele,:]], dims=2)
    ce_p_gs = sum(gs_pe.Ni .* ce_p[element_pe[gs_pe.ele,:]], dims=2)
    ce_sp_gs = sum(gs_sp.Ni .* ce_sp[element_sp[gs_sp.ele, :]], dims = 2)
    u_n_gs = param.NE.U(csn_surf_gs)
    u_p_gs = param.PE.U(csp_surf_gs)
    eta_p_gs = phis_p[1] .- phie_p_gs - u_p_gs
    eta_n_gs = phis_n[1] .- phie_n_gs - u_n_gs
    j0_n_gs =  param.NE.k * Arrhenius(param.NE.Eac_k, T) .* abs.(csn_surf_gs .* (1.0 .- csn_surf_gs) .* ce_n_gs) .^ 0.5
    j0_p_gs =  param.PE.k * Arrhenius(param.PE.Eac_k, T) .* abs.(csp_surf_gs .* (1.0 .- csp_surf_gs) .* ce_p_gs) .^ 0.5
    j_n_gs = j0_n_gs .* sinh.(0.5 * eta_n_gs / T) * 2.0
    j_p_gs = j0_p_gs .* sinh.(0.5 * eta_p_gs / T) * 2.0

    variables["negative electrode interfacial current density"] = j_n
    variables["positive electrode interfacial current density"] = j_p
    variables["negative electrode exchange current density"] = j0_n
    variables["positive electrode exchange current density"] = j0_p
    variables["negative electrode open circuit potential"] = u_n
    variables["positive electrode open circuit potential"] = u_p
    variables["negative electrode overpotential"] = eta_n
    variables["positive electrode overpotential"] = eta_p    
    variables["negative electrode interfacial current at Gauss point"] = j_n_gs
    variables["positive electrode interfacial current at Gauss point"] = j_p_gs
    variables["negative electrode exchange current density at Gauss point"] = j0_n_gs
    variables["positive electrode exchange current density at Gauss point"] = j0_p_gs
    variables["negative electrode open circuit potential at Gauss point"] = u_n_gs
    variables["positive electrode open circuit potential at Gauss point"] = u_p_gs
    variables["negative electrode overpotential at Gauss point"] = eta_n_gs
    variables["positive electrode overpotential at Gauss point"] = eta_p_gs
    variables["electrolyte lithium concentration at negative electrode Gauss point"] = ce_n_gs
    variables["electrolyte lithium concentration at positive electrode Gauss point"] = ce_p_gs
    variables["electrolyte lithium concentration at separator Gauss point"] = ce_sp_gs
    variables["time"] = t
    variables["temperature"] = T
    variables["cell voltage"] = phis_p[end] - phis_n[1]
    variables["cell current"] =case.opt.Current(t * case.param.scale.t0) / case.param_dim.cell.I1C
    return variables
end
