function SPM(case::Case)

    if case.opt.jacobi_K == "constant"
        K = []
        M = []
    else
        param = case.param
        mesh1 = case.mesh["negative particle"]
        mesh2 = case.mesh["positive particle"]
        M1, K1 = ElectrodeDiffusion(param.NE, mesh1, mesh1.nlen, case.opt)
        M2, K2 = ElectrodeDiffusion(param.PE, mesh2, mesh2.nlen, case.opt)   
        M1 .*= param.scale.ts_n
        M2 .*= param.scale.ts_p 

        K = blockdiag(K1, K2)
        M = blockdiag(M1, M2)
    end
    t=0.1
    F = SPM_BC(case, t)

    return M, K, F
end


function SPM_BC(case::Case, t::Float64)
    param = case.param
    I = case.opt.Current(t) / case.param_dim.cell.area / param.scale.I_typ

    flux1 = zeros(Float64, case.mesh["negative particle"].nlen, 1)
    flux1[end] = I / param.NE.as / param.NE.thickness * param.NE.rs^2

    flux2 = zeros(Float64, case.mesh["positive particle"].nlen, 1)
    flux2[end] = - I / param.PE.as / param.PE.thickness * param.PE.rs^2

    F = [flux1; flux2]
    return F
end

function SPM_update!(case::Case, variables::Dict{String, Matrix{Float64}}, v::Int64, yt::Array{Float64}, t::Float64, T::Float64=1.0)
    param = case.param
    I = case.opt.Current(t) / case.param_dim.cell.area / param.scale.I_typ

    cn_surf = yt[case.index["csn"][end], 1]
    cp_surf = yt[case.index["csp"][end], 1]
    u_n = param.NE.U(cn_surf)
    u_p = param.PE.U(cp_surf)
    j0_n =  param.NE.k * Arrhenius(param.NE.Eac_k, T) * sqrt(cn_surf * param.EL.ce0 * (1 - cn_surf))
    j0_p =  param.PE.k * Arrhenius(param.PE.Eac_k, T) * sqrt(cp_surf * param.EL.ce0 * (1 - cp_surf))
    eta_n = 2 * T * asinh(-I / param.NE.as / param.NE.thickness /2 / j0_n)
    eta_p = 2 * T * asinh(I / param.PE.as / param.PE.thickness / 2 / j0_p)
    V_cell = u_p - u_n + eta_p - eta_n  # overpotentials seem not right
    variables["negative particle surface lithium concentration"][1,v] = cn_surf
    variables["positive particle surface lithium concentration"][1,v] = cp_surf
    variables["cell voltage"][1,v] = V_cell 
    variables["negative electrode exchange current density"][1,v] = j0_n
    variables["positive electrode exchange current density"][1,v] = j0_p
    variables["negative electrode overpotential"][1,v] = eta_n
    variables["positive electrode overpotential"][1,v] = eta_p
    variables["time"][1,v] = t
end