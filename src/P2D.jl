function P2D(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)

    variables = P2D_variables(case, yt, t)
    F = P2D_BC(case, variables)
    if jacobi == "constant" # no need to update M and K
        M = []
        K = []
    else
        param = case.param
        mesh1 = case.mesh["negative particle"]
        mesh2 = case.mesh["positive particle"]
        mesh3 = case.mesh["electrolyte"]
        mesh4 = case.mesh["negative electrode"]
        mesh5 = case.mesh["positive electrode"]
        M1, K1 = ElectrodeDiffusion(param.NE, mesh1, mesh1.nlen)
        M2, K2 = ElectrodeDiffusion(param.PE, mesh2, mesh2.nlen) 
        M3, K3 = ElectrolyteDiffusion(param, mesh3, mesh3.nlen, variables)   
        M4, K4 = ElectrodePotential(param.NE, mesh4, mesh4.nlen)
        M5, K5 = ElectrodePotential(param.PE, mesh5, mesh5.nlen) 
        M6, K6 = ElectrolytePotential(param, mesh3, mesh3.nlen, variables) 
        K4[1,2:end] .= 0.0
        K4[1,1] = 1.0
        F[size(K1,1) + size(K2,1) + size(K3,1) + 1] = 0.0
        M1 = M1 .* param.scale.ts_n
        M2 = M2 .* param.scale.ts_p
        M3 = M3 .* param.scale.te 

        # # the following part seems not to affect the result, need to recheck later
        # # add interface boundary condition 
        # xloc = [1.0, 0, 1.0, 0]
        # v = [case.opt.Nn, case.opt.Nn + 1, case.opt.Nn + case.opt.Ns, case.opt.Nn + case.opt.Ns + 1]
        # _, dNidx = ShapeFunction1D(mesh3.element, mesh3.type, mesh3.node, xloc, v)  
        # v_ns = case.mesh["negative electrode"].nlen
        # v_sp = case.mesh["negative electrode"].nlen + case.mesh["separator"].nlen - 1
        # M3[v_ns, :] .= 0
        # M3[v_sp, :] .= 0
        # K3[v_ns, :] .= 0
        # K3[v_sp, :] .= 0
        # ce = param.EL.ce0 # need to modify later
        # T = 298.0
        # De_eff =  param.EL.De(ce) * Arrhenius(param.EL.Eac_D, T)
        # K3[v_ns, mesh3.element[v[1],:]] .+= dNidx[1,:] .* De_eff
        # K3[v_ns, mesh3.element[v[2],:]] .+= - dNidx[2,:] .* De_eff
        # K3[v_sp, mesh3.element[v[3],:]] .+= dNidx[3,:] .* De_eff
        # K3[v_sp, mesh3.element[v[4],:]] .+= - dNidx[4,:] .* De_eff
        # F[mesh1.nlen + mesh2.nlen + v_ns] = 0
        # F[mesh1.nlen + mesh2.nlen + v_sp] = 0


        M = blockdiag(M1, M2, M3, M4, M5, M6)
        K = blockdiag(K1, K2, K3, K4, K5, K6)
    end

    return M, K, F, variables
end

function P2D_BC(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    t = variables["time"]
    I = case.opt.Current(t) / case.param_dim.cell.area / param.scale.I_typ

    # electrode particle diffusion source term
    j_n = variables["negative electrode interfacial current density"]
    j_p = variables["positive electrode interfacial current density"]
    flux_np = zeros(Float64, case.mesh["negative particle"].nlen, 1)
    v_np_surf = case.index["negative particle surface lithium concentration"]
    flux_np[v_np_surf] = j_n * param.NE.rs^2
    flux_pp = zeros(Float64, case.mesh["positive particle"].nlen, 1)
    v_pp_surf = case.index["positive particle surface lithium concentration"] .- case.mesh["negative particle"].nlen
    flux_pp[v_pp_surf] = j_p * param.PE.rs^2

    # electrolyte diffusion source term
    mesh_el = case.mesh["electrolyte"]
    n_sp_gs = case.opt.Ns * mesh_el.gs.order
    j_n_gs = variables["negative electrode interfacial current at Gauss point"]
    j_p_gs = variables["positive electrode interfacial current at Gauss point"]
    j_el_gs = [param.NE.as .* j_n_gs; zeros(n_sp_gs, 1); param.PE.as .* j_p_gs] 
    coeff_el = mesh_el.gs.weight .* mesh_el.gs.detJ .* j_el_gs * (1 - param.EL.tplus)
    Vi_el = mesh_el.element[mesh_el.gs.ele,:]
    flux_el = Assemble1D(Vi_el, mesh_el.gs.Ni, coeff_el, mesh_el.nlen)

    # negative electrode charge source term
    mesh_ne = case.mesh["negative electrode"]
    coeff_ne = mesh_ne.gs.weight .* mesh_ne.gs.detJ * param.NE.as .* j_n_gs
    Vi_ne = mesh_ne.element[mesh_ne.gs.ele,:]
    flux_ne = Assemble1D(Vi_ne, mesh_ne.gs.Ni, coeff_ne, mesh_ne.nlen)
    flux_ne[1] -= I

    # positive electrode charge source term
    mesh_pe = case.mesh["positive electrode"]
    coeff_pe = mesh_pe.gs.weight .* mesh_pe.gs.detJ * param.PE.as .* j_p_gs
    Vi_pe = mesh_pe.element[mesh_pe.gs.ele,:]
    flux_pe = Assemble1D(Vi_pe, mesh_pe.gs.Ni, coeff_pe, mesh_pe.nlen)
    flux_pe[end] -= I

    # electrolyte charge source term
    coeff_elc = mesh_el.gs.weight .* mesh_el.gs.detJ .* j_el_gs
    flux_elc = Assemble1D(Vi_el, mesh_el.gs.Ni, coeff_elc, mesh_el.nlen)

    ce = variables["electrolyte lithium concentration"]
    ce_gs = sum(mesh_el.gs.Ni .* ce[mesh_el.element[mesh_el.gs.ele,:]], dims=2)
    dcedx_gs = sum(mesh_el.gs.dNidx .* ce[mesh_el.element[mesh_el.gs.ele,:]], dims=2)
    T = variables["temperature"]
    kappa_D_eff = param.EL.kappa(ce_gs) .* mesh_el.gs.weight .* mesh_el.gs.detJ
    tau_ne = param.NE.eps ^ param.NE.brugg * ones(size(j_n_gs))
    tau_pe = param.PE.eps ^ param.PE.brugg * ones(size(j_p_gs))
    tau_sp = param.SP.eps ^ param.SP.brugg * ones(n_sp_gs, 1)
    tau_el = [tau_ne; tau_sp; tau_pe]
    kappa_D_eff .*= 2 * T * (1 - param.EL.tplus) .* param.EL.dlnf_dlnc(ce) .* tau_el
    kappa_D_eff .*= 1 ./ ce_gs .* dcedx_gs
    flux_elc .+= Assemble1D(Vi_el, mesh_el.gs.dNidx, kappa_D_eff, mesh_el.nlen)
    # assemble all fluxes
    flux = [flux_np; flux_pp; flux_el; flux_ne; flux_pe; flux_elc]
    return flux
end

function P2D_variables(case::Case, yt::Array{Float64}, t::Float64)
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
    j0_n =  param.NE.k * Arrhenius(param.NE.Eac_k, T) .* (csn_surf .* (1.0 .- csn_surf) .* ce_n) .^ 0.5
    j0_p =  param.PE.k * Arrhenius(param.PE.Eac_k, T) .* (csp_surf .* (1.0 .- csp_surf) .* ce_p) .^ 0.5
    j_n = j0_n .* sinh.(0.5 .* eta_n) * 2.0
    j_p = j0_p .* sinh.(0.5 .* eta_p) * 2.0

    csn_surf_gs = sum(gs_ne.Ni .* csn_surf[element_ne[gs_ne.ele,:]], dims=2)
    csp_surf_gs = sum(gs_pe.Ni .* csp_surf[element_pe[gs_pe.ele,:]], dims=2)
    phis_n_gs = sum(gs_ne.Ni .* phis_n[element_ne[gs_ne.ele,:]], dims=2)
    phis_p_gs = sum(gs_pe.Ni .* phis_p[element_pe[gs_pe.ele,:]], dims=2)
    phie_n_gs = sum(gs_ne.Ni .* phie_n[element_ne[gs_ne.ele,:]], dims=2)
    phie_p_gs = sum(gs_pe.Ni .* phie_p[element_pe[gs_pe.ele,:]], dims=2)
    ce_n_gs = sum(gs_ne.Ni .* ce_n[element_ne[gs_ne.ele,:]], dims=2)
    ce_p_gs = sum(gs_pe.Ni .* ce_p[element_pe[gs_pe.ele,:]], dims=2)
    ce_sp_gs = sum(gs_sp.Ni .* ce_sp[element_sp[gs_sp.ele, :]], dims = 2)
    u_n_gs = param.NE.U(csn_surf_gs)
    u_p_gs = param.PE.U(csp_surf_gs)
    eta_p_gs = phis_p_gs - phie_p_gs - u_p_gs
    eta_n_gs = phis_n_gs - phie_n_gs - u_n_gs
    j0_n_gs =  param.NE.k * Arrhenius(param.NE.Eac_k, T) .* (csn_surf_gs .* (1.0 .- csn_surf_gs) .* ce_n_gs) .^ 0.5
    j0_p_gs =  param.PE.k * Arrhenius(param.PE.Eac_k, T) .* (csp_surf_gs .* (1.0 .- csp_surf_gs) .* ce_p_gs) .^ 0.5
    j_n_gs = j0_n_gs .* sinh.(0.5 * eta_n_gs) * 2.0
    j_p_gs = j0_p_gs .* sinh.(0.5 * eta_p_gs) * 2.0

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
    return variables
end