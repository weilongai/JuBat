function SPM(case::Case, y0::Array{Float64}, t::Float64, variables::Dict{String, Any})
    # negative elecrode particle diffusion 
    param = case.param
    mesh1 = case.mesh["negative particle"]
    I = case.opt.Current(t) / case.param_dim.cell.Total_surface / param.scale.I_typ
    flux1 = zeros(Float64, mesh1.nlen, 1)
    flux1[end] = - I / param.NE.as / param.NE.thickness 

    M1, K1, F1 = ElectrodeDiffusion(param.NE, mesh1, flux1, mesh1.nlen)
    M1 = M1 .* param.scale.ts_n

    # positive elecrode particle diffusion
    mesh2 = case.mesh["positive particle"]
    flux2 = zeros(Float64, mesh2.nlen, 1)
    flux2[end] = - I / param.PE.as / param.PE.thickness 
    M2, K2, F2 = ElectrodeDiffusion(param.PE, mesh2, flux2, mesh2.nlen)
    M2 = M2 .* param.scale.ts_p 
    M = blockdiag(M1, M2)
    K = blockdiag(K1, K2)
    F = [F1; F2]

    vt = variables["time step"]
    cn_surf = y0[mesh1.nlen, 1]
    cp_surf = y0[mesh1.nlen + mesh2.nlen, 1]
    u_n = param.NE.OCP(cn_surf)
    u_p = param.PE.OCP(cp_surf)
    j0_n = param.NE.k * sqrt(cn_surf * param.EL.ce0 * (1 - cn_surf))
    j0_p = param.PE.k * sqrt(cp_surf * param.EL.ce0 * (1 - cp_surf))
    V_cell = u_p - u_n - 2 * asinh(I / param.PE.as / param.PE.thickness / j0_p) - 2 * asinh(I / param.NE.as / param.NE.thickness / j0_n)
    variables["negative particle surface lithium concentration"][1,vt] = cn_surf
    variables["positive particle surface lithium concentration"][1,vt] = cp_surf
    variables["cell voltage"][1,vt] = V_cell 
    variables["negative electrode exchange current density"][1,vt] = j0_n
    variables["positive electrode exchange current density"][1,vt] = j0_p
    # now need to solve Ma=Kc+f
    return M, K, F, variables
end