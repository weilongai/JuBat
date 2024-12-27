function SPM(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    variables = SPM_variables(case, yt, t)
    if case.opt.mechanicalmodel == "full"
        variables = Mechanicaloutput(case,variables)
        theta_Mn = variables["negative particle stress coupling diffusion coefficient"][1]
        theta_Mp = variables["positive particle stress coupling diffusion coefficient"][1]
        else
        theta_Mn = 0.0
        theta_Mp = 0.0
    end
    csn_gs = variables["negative particle concentration at gauss point"]
    csp_gs = variables["positive particle concentration at gauss point"]
    if jacobi == "constant" && param.NE.M_d != [] # no need to update M and K
        M_np = param.NE.M_d
        K_np = param.NE.K_d
        M_pp = param.PE.M_d
        K_pp = param.PE.K_d
    else
        param = case.param
        mesh_np = case.mesh["negative particle"]
        mesh_pp = case.mesh["positive particle"]
        M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, theta_Mn)
        M_pp, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, theta_Mp)   
        M_np .*= param.scale.ts_n / param_dim.scale.t0
        M_pp .*= param.scale.ts_p / param_dim.scale.t0
    end
    K = blockdiag(K_np, K_pp)
    M = blockdiag(M_np, M_pp)
    F = SPM_BC(case, variables)
    return M, K, F, variables
end


function SPM_BC(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    j_n = variables["negative electrode interfacial current density"]
    j_p = variables["positive electrode interfacial current density"] 

    flux_np = zeros(Float64, case.mesh["negative particle"].nlen, 1)
    flux_np[end] = - j_n * param.NE.rs^2

    flux_pp = zeros(Float64, case.mesh["positive particle"].nlen, 1)
    flux_pp[end] = - j_p * param.PE.rs^2

    F = [flux_np; flux_pp]
    return F
end

function SPM_variables(case::Case, yt::Array{Float64}, t::Float64)
    param = case.param
    I_app =case.opt.Current(t * case.param.scale.t0) / param.scale.I_typ
    j_n = I_app / param.NE.as / param.NE.thickness
    j_p = - I_app / param.PE.as / param.PE.thickness
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
    cn_surf = variables["negative particle surface lithium concentration"]
    cp_surf = variables["positive particle surface lithium concentration"]
    u_n = param.NE.U(cn_surf) .+ (T .- case.param.cell.T0) * param.NE.dUdT(cn_surf)
    u_p = param.PE.U(cp_surf) .+ (T .- case.param.cell.T0) * param.PE.dUdT(cp_surf)
    j0_n =  param.NE.k * Arrhenius(param.NE.Eac_k, T) * sqrt.(cn_surf .* param.EL.ce0 .* abs.(1.0 .- cn_surf))
    j0_p =  param.PE.k * Arrhenius(param.PE.Eac_k, T) * sqrt.(cp_surf .* param.EL.ce0 .* abs.(1.0 .- cp_surf))
    eta_n = 2 * T .* asinh.(j_n / 2.0 ./ j0_n)
    eta_p = 2 * T .* asinh.(j_p / 2.0 ./ j0_p)
    V_cell = u_p - u_n + eta_p - eta_n 
    variables["cell voltage"] = V_cell[1]
    variables["negative electrode exchange current density"] = j0_n
    variables["positive electrode exchange current density"] = j0_p
    variables["negative electrode interfacial current density"] = j_n
    variables["positive electrode interfacial current density"] = j_p
    variables["negative electrode overpotential"] = eta_n
    variables["positive electrode overpotential"] = eta_p
    variables["negative electrode open circuit potential"] = u_n
    variables["positive electrode open circuit potential"] = u_p
    variables["time"] = t
    variables["temperature"] = T
    variables["cell current"] =case.opt.Current(t * case.param.scale.t0) / case.param_dim.cell.I1C
    return variables
end