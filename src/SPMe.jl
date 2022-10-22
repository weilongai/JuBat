function SPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    variables = SPMe_variables(case, yt, t)
    param = case.param
    if jacobi == "constant" && param.NE.M_d != [] # no need to update M and K
        M_np = param.NE.M_d
        K_np = param.NE.K_d
        M_pp = param.PE.M_d
        K_pp = param.PE.K_d
    else
        mesh_np = case.mesh["negative particle"]
        mesh_pp = case.mesh["positive particle"]
        M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen)
        M_pp, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen) 
        M_np = M_np .* param.scale.ts_n
        M_pp = M_pp .* param.scale.ts_p
    end
    mesh_el = case.mesh["electrolyte"]        
    M_el, K_el = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, variables)   
    M_el = M_el .* param.scale.te 
    F = SPMe_BC(case, variables)

    # # the following part seems not to affect the result, need to recheck later
    # # add interface boundary condition 
    # xloc = [1.0, 0, 1.0, 0]
    # v = [case.opt.Nn, case.opt.Nn + 1, case.opt.Nn + case.opt.Ns, case.opt.Nn + case.opt.Ns + 1]
    # _, dNidx = ShapeFunction1D(mesh_el.element, mesh_el.type, mesh_el.node, xloc, v)  
    # v_ns = case.mesh["negative electrode"].nlen
    # v_sp = case.mesh["negative electrode"].nlen + case.mesh["separator"].nlen - 1
    # M_el[v_ns, :] .= 0
    # M_el[v_sp, :] .= 0
    # K_el[v_ns, :] .= 0
    # K_el[v_sp, :] .= 0
    # ce = param.EL.ce0 # need to modify later
    # T = variables["temperature"]
    # De_eff =  param.EL.De(ce) * Arrhenius(param.EL.Eac_D, T)
    # K_el[v_ns, mesh_el.element[v[1],:]] .+= dNidx[1,:] .* De_eff * param.NE.eps ^ param.NE.brugg
    # K_el[v_ns, mesh_el.element[v[2],:]] .+= - dNidx[2,:] .* De_eff * param.SP.eps ^ param.SP.brugg
    # K_el[v_sp, mesh_el.element[v[3],:]] .+= dNidx[3,:] .* De_eff * param.SP.eps ^ param.SP.brugg
    # K_el[v_sp, mesh_el.element[v[4],:]] .+= - dNidx[4,:] .* De_eff * param.PE.eps ^ param.PE.brugg
    # F[mesh_np.nlen + mesh_pp.nlen + v_ns] = 0
    # F[mesh_np.nlen + mesh_pp.nlen + v_sp] = 0
    
    M = blockdiag(M_np, M_pp, M_el)
    K = blockdiag(K_np, K_pp, K_el)

    return M, K, F, variables
end


function SPMe_BC(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    j_n = variables["negative electrode interfacial current density"]
    j_p = variables["positive electrode interfacial current density"]
    flux_np = zeros(Float64, case.mesh["negative particle"].nlen, 1)
    flux_np[end] = - j_n * param.NE.rs^2
    flux_pp = zeros(Float64, case.mesh["positive particle"].nlen, 1)
    flux_pp[end] = - j_p * param.PE.rs^2

    # electrolyte source term
    mesh_el = case.mesh["electrolyte"]
    coeff = mesh_el.gs.weight .* mesh_el.gs.detJ
    v_ne = collect(1:case.opt.Nn * mesh_el.gs.order)
    v_sp = case.opt.Nn * mesh_el.gs.order .+ collect(1:case.opt.Ns * mesh_el.gs.order)
    v_pe = (case.opt.Nn + case.opt.Ns) * mesh_el.gs.order .+ collect(1:case.opt.Np * mesh_el.gs.order)
    coeff[v_ne] .*= (1 - param.EL.tplus) * param.NE.as .* j_n
    coeff[v_sp] .*= 0
    coeff[v_pe] .*= (1 - param.EL.tplus) * param.PE.as .* j_p

    Vi = mesh_el.element[mesh_el.gs.ele,:]
    flux_el = Assemble1D(Vi, mesh_el.gs.Ni, coeff, mesh_el.nlen)
    flux = [flux_np; flux_pp; flux_el]
    return flux
end

function SPMe_variables(case::Case, yt::Array{Float64}, t::Float64)
    param = case.param
    variables = StandardVariables(case, 1)
    I_app = case.opt.Current(t) / case.param_dim.cell.area / param.scale.I_typ
    j_n = - I_app / param.NE.as / param.NE.thickness
    j_p = I_app / param.PE.as / param.PE.thickness
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    mesh_sp = case.mesh["separator"]
    var_list = collect(keys(case.index))
    for i in var_list
        variables[i] = yt[case.index[i]]
    end
    if "temperature" in var_list
        T = yt[case.index["temperature"]]
    else
        T = case.param.cell.T0
    end
    cn_surf = variables["negative particle surface lithium concentration"]
    cp_surf = variables["positive particle surface lithium concentration"]
    ce_n = variables["electrolyte lithium concentration in negative electrode"]
    ce_p = variables["electrolyte lithium concentration in positive electrode"]
    ce_sp = variables["electrolyte lithium concentration in separator"]

    ce_n_gs = sum(mesh_ne.gs.Ni .* ce_n[mesh_ne.element[mesh_ne.gs.ele, :]], dims = 2)
    ce_p_gs = sum(mesh_pe.gs.Ni .* ce_p[mesh_pe.element[mesh_pe.gs.ele, :]], dims = 2)
    ce_sp_gs = sum(mesh_sp.gs.Ni .* ce_sp[mesh_sp.element[mesh_sp.gs.ele, :]], dims = 2)

    j0_n_gs =  param.NE.k * Arrhenius(param.NE.Eac_k, T) .* (cn_surf .* abs.(1.0 .- cn_surf) .* ce_n_gs) .^ 0.5
    j0_p_gs =  param.PE.k * Arrhenius(param.PE.Eac_k, T) .* (cp_surf .* abs(1.0 .- cp_surf) .* ce_p_gs) .^ 0.5
    j0_n_av = IntV(j0_n_gs, mesh_ne) / param.NE.thickness
    j0_p_av = IntV(j0_p_gs, mesh_pe) / param.PE.thickness
    eta_n = 2.0 * T * asinh.(j_n / 2.0 / j0_n_av)
    eta_p = 2.0 * T * asinh.(j_p / 2.0 / j0_p_av)

    # kappa_ne = IntV(param.EL.kappa(ce_n_gs), mesh_ne) / param.NE.thickness * param.NE.eps ^ param.NE.brugg
    # kappa_pe = IntV(param.EL.kappa(ce_p_gs), mesh_pe) / param.PE.thickness * param.PE.eps ^ param.PE.brugg
    # kappa_sp = IntV(param.EL.kappa(ce_sp_gs), mesh_sp) / param.SP.thickness * param.SP.eps ^ param.SP.brugg
    # R_EL = param.NE.thickness / kappa_ne  + 2 * param.SP.thickness / kappa_sp + param.PE.thickness / kappa_pe
    # dphi_e = 2.0 * T * (1 - param.EL.tplus) * log(ce_p[end]/ce_n[1]) + I_app / 2.0 * R_EL # it seems to be minus here

    ## below expression by Ferran et al. progress in energy 2022 gives better accuracy 
    kappa_ne = param.EL.kappa(param.EL.ce0) * param.NE.eps ^ param.NE.brugg
    kappa_pe = param.EL.kappa(param.EL.ce0) * param.PE.eps ^ param.PE.brugg
    kappa_sp = param.EL.kappa(param.EL.ce0) * param.SP.eps ^ param.SP.brugg
    R_EL = param.NE.thickness / kappa_ne / 3.0  + param.SP.thickness / kappa_sp + param.PE.thickness / kappa_pe / 3.0
    dphi_e = 2.0 * T * (1 - param.EL.tplus) * log(ce_p[end]/ce_n[1]) + I_app * R_EL 


    u_n = param.NE.U(cn_surf)
    u_p = param.PE.U(cp_surf)   
    V_cell = u_p - u_n + eta_p - eta_n + dphi_e # + eta_S
    variables["negative particle surface lithium concentration"] = cn_surf
    variables["positive particle surface lithium concentration"] = cp_surf
    variables["cell voltage"] = V_cell 
    variables["negative electrode exchange current density"] = j0_n_av
    variables["positive electrode exchange current density"] = j0_p_av
    variables["negative electrode interfacial current density"] = j_n
    variables["positive electrode interfacial current density"] = j_p
    variables["negative electrode overpotential"] = eta_n
    variables["positive electrode overpotential"] = eta_p
    variables["negative electrode open circuit potential"] = u_n
    variables["positive electrode open circuit potential"] = u_p
    variables["electrolyte lithium concentration at negative electrode Gauss point"] = ce_n_gs
    variables["electrolyte lithium concentration at positive electrode Gauss point"] = ce_p_gs
    variables["electrolyte lithium concentration at separator Gauss point"] = ce_sp_gs
    variables["time"] = t
    variables["temperature"] = T
    variables["cell current"] = case.opt.Current(t) / case.param_dim.cell.I1C
    return variables
end