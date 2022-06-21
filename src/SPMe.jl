function SPMe(case::Case, opt::String)
    param = case.param
    mesh1 = case.mesh["negative particle"]
    mesh2 = case.mesh["positive particle"]
    mesh3 = case.mesh["electrolyte"]

    K1 = ElectrodeDiffusion(param.NE, mesh1, mesh1.nlen, opt)
    K2 = ElectrodeDiffusion(param.PE, mesh2, mesh2.nlen, opt) 
    K3 = ElectrolyteDiffusion(param, mesh3, mesh3.nlen, opt)   
    if opt == "M"
        K1 = K1 .* param.scale.ts_n
        K2 = K2 .* param.scale.ts_p
        K3 = K3 .* param.scale.te 
    # elseif opt == "K"
    #     # add interface boundary condition
    #     xloc = [1.0, 0, 1.0, 0]
    #     v = [case.opt.Nn, case.opt.Nn, case.opt.Nn + case.opt.Ns, case.opt.Nn + case.opt.Ns + 1]
    #     _, dNidx = ShapeFunction1D(mesh3.element, mesh3.type, mesh3.node, xloc, v)  
    #     v_ns = case.mesh["negative electrode"].nlen
    #     v_sp = case.mesh["negative electrode"].nlen + case.mesh["separator"].nlen - 1
    #     K3[v_ns, :] .= 0
    #     K3[v_sp, :] .= 0
    #     ce = param.EL.ce0 # need to modify later
    #     T = 298.0
    #     De_eff =  param.EL.De(ce) * Arrhenius(param.EL.Eac_D, T)
    #     K3[v_ns, mesh3.element[case.opt.Nn,:]] .+= dNidx[1,:] .* De_eff
    #     K3[v_ns, mesh3.element[case.opt.Nn + 1,:]] .+= - dNidx[2,:] .* De_eff
    #     K3[v_sp, mesh3.element[case.opt.Nn + case.opt.Ns, :]] .+= dNidx[3,:] .* De_eff
    #     K3[v_sp, mesh3.element[case.opt.Nn + case.opt.Ns + 1,:]] .+= - dNidx[4,:] .* De_eff
    end
    K = blockdiag(K1, K2, K3)
    return K
end


function SPMe_BC(case::Case, t::Float64)
    param = case.param
    I = case.opt.Current(t) / case.param_dim.cell.area / param.scale.I_typ
    j_n = I / param.NE.as / param.NE.thickness 
    j_p = - I / param.PE.as / param.PE.thickness
    flux_np = zeros(Float64, case.mesh["negative particle"].nlen, 1)
    flux_np[end] = j_n * param.NE.Rs^2
    flux_pp = zeros(Float64, case.mesh["positive particle"].nlen, 1)
    flux_pp[end] = j_p * param.PE.Rs^2

    # electrolyte source term
    EL_mesh = case.mesh["electrolyte"]
    coeff = EL_mesh.gs.weight .* EL_mesh.gs.detJ
    v_ne = collect(1:case.opt.Nn * EL_mesh.gs.order)
    v_sp = case.opt.Nn * EL_mesh.gs.order .+ collect(1:case.opt.Ns * EL_mesh.gs.order)
    v_pe = (case.opt.Nn + case.opt.Ns) * EL_mesh.gs.order .+ collect(1:case.opt.Np * EL_mesh.gs.order)
    coeff[v_ne] .*= (1 - param.EL.tplus) * param.NE.as * j_n
    coeff[v_sp] .*= 0
    coeff[v_pe] .*= (1 - param.EL.tplus) * param.PE.as * j_p

    Vi = EL_mesh.element[EL_mesh.gs.ele,:]
    flux_el = Assemble1D(Vi, EL_mesh.gs.Ni, coeff, EL_mesh.nlen)

    # flux_el[case.mesh["negative electrode"].nlen] = 0
    # flux_el[case.mesh["negative electrode"].nlen + case.mesh["separator"].nlen - 1] = 0
    flux = [flux_np; flux_pp; flux_el]
    return flux
end

function SPMe_update!(case::Case, variables::Dict{String, Matrix{Float64}}, v::Int64, yt::Array{Float64}, t::Float64, T::Float64=1.)
    param = case.param
    I = case.opt.Current(t) / case.param_dim.cell.area / param.scale.I_typ
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    mesh_sp = case.mesh["separator"]
    cn_surf = yt[case.index["csn"][end], 1]
    cp_surf = yt[case.index["csp"][end], 1]


    ce = yt[case.index["ce"],1]
    ce_n_gs = sum(mesh_ne.gs.Ni .* ce[mesh_ne.element[mesh_ne.gs.ele, :]], dims = 2)
    ce_p_gs = sum(mesh_pe.gs.Ni .* ce[mesh_pe.element[mesh_pe.gs.ele, :]], dims = 2)
    ce_sp_gs = sum(mesh_sp.gs.Ni .* ce[mesh_sp.element[mesh_sp.gs.ele, :]], dims = 2)

    j0_n_gs =  param.NE.k * Arrhenius(param.NE.Eac_k, T) .* (cn_surf * (1 - cn_surf) .* ce_n_gs) .^ 0.5
    j0_p_gs =  param.PE.k * Arrhenius(param.PE.Eac_k, T) .* (cp_surf * (1 - cp_surf) .* ce_p_gs) .^ 0.5
    j0_n_av = IntV(j0_n_gs, mesh_ne) / param.NE.thickness
    j0_p_av = IntV(j0_p_gs, mesh_pe) / param.PE.thickness
    eta_n = 2 * T * asinh(-I / param.NE.as / param.NE.thickness /2 / j0_n_av)
    eta_p = 2 * T * asinh(I / param.PE.as / param.PE.thickness / 2 / j0_p_av)

    kappa_ne = IntV(param.EL.kappa(ce_n_gs), mesh_ne) / param.NE.thickness
    kappa_pe = IntV(param.EL.kappa(ce_p_gs), mesh_pe) / param.PE.thickness
    kappa_sp = IntV(param.EL.kappa(ce_sp_gs), mesh_sp) / param.SP.thickness
    R_EL = param.NE.thickness / kappa_ne  + 2* param.SP.thickness / kappa_sp + param.PE.thickness / kappa_pe
    eta_EL = 2 * T * param.EL.tplus * log(abs(ce[end]/ce[1])) - I / 2 * R_EL

    u_n = param.NE.U(cn_surf)
    u_p = param.PE.U(cp_surf)   
    V_cell = u_p - u_n + eta_p - eta_n + eta_EL  # overpotentials seem not right
    variables["negative particle surface lithium concentration"][v] = cn_surf
    variables["positive particle surface lithium concentration"][v] = cp_surf
    variables["cell voltage"][v] = V_cell 
    variables["negative electrode exchange current density"][v] = j0_n_av
    variables["positive electrode exchange current density"][v] = j0_p_av
    variables["negative electrode overpotential"][v] = eta_n
    variables["positive electrode overpotential"][v] = eta_p
    variables["time"][v] = t
    variables["electrolyte lithium concentration"][:,v] = ce
end