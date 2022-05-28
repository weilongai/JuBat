function SPM(case::Case, opt::String)
    param = case.param
    mesh1 = case.mesh["negative particle"]
    mesh2 = case.mesh["positive particle"]

    K1 = ElectrodeDiffusion(param.NE, mesh1, mesh1.nlen, opt)
    K2 = ElectrodeDiffusion(param.PE, mesh2, mesh2.nlen, opt)   
    if opt == "M"
        K1 = K1 .* param.scale.ts_n
        K2 = K2 .* param.scale.ts_p 
    end
    K = blockdiag(K1, K2)
    return K
end


function SPM_BC(case::Case, t::Float64)
    param = case.param
    I = case.opt.Current(t) / case.param_dim.cell.Total_surface / param.scale.I_typ

    flux1 = zeros(Float64, case.mesh["negative particle"].nlen, 1)
    flux1[end] = I / param.NE.as / param.NE.thickness 

    flux2 = zeros(Float64, case.mesh["positive particle"].nlen, 1)
    flux2[end] = - I / param.PE.as / param.PE.thickness 

    F = [flux1 * param.NE.Rs^2; flux2 * param.PE.Rs^2]
    return F
end

function SPM_update!(case::Case, variables::Dict{String, Matrix{Float64}}, v::Integer, yt::Array{Float64}, t::Float64)
    param = case.param
    I = case.opt.Current(t) / case.param_dim.cell.Total_surface / param.scale.I_typ

    cn_surf = yt[case.index["csn"][end], 1]
    cp_surf = yt[case.index["csp"][end], 1]
    u_n = param.NE.OCP(cn_surf)
    u_p = param.PE.OCP(cp_surf)
    j0_n =  param.NE.k * sqrt(cn_surf * param.EL.ce0 * (1 - cn_surf))
    j0_p =  param.PE.k * sqrt(cp_surf * param.EL.ce0 * (1 - cp_surf))
    V_cell = u_p - u_n + 2 * asinh(I / param.PE.as / param.PE.thickness / 2 / j0_p) + 2 * asinh(I / param.NE.as / param.NE.thickness /2 / j0_n) # overpotentials seem not right
    variables["negative particle surface lithium concentration"][1,v] = cn_surf
    variables["positive particle surface lithium concentration"][1,v] = cp_surf
    variables["cell voltage"][1,v] = V_cell 
    variables["negative electrode exchange current density"][1,v] = j0_n
    variables["positive electrode exchange current density"][1,v] = j0_p
    variables["time"][1,v] = t
end