function SPM(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    variables = SPM_variables(case, yt, t)
    if jacobi == "constant" || jacobi == "constant MK" 
        K = []
        M = []
    else
        param = case.param
        mesh1 = case.mesh["negative particle"]
        mesh2 = case.mesh["positive particle"]
        M1, K1 = ElectrodeDiffusion(param.NE, mesh1, mesh1.nlen)
        M2, K2 = ElectrodeDiffusion(param.PE, mesh2, mesh2.nlen)   
        M1 .*= param.scale.ts_n
        M2 .*= param.scale.ts_p 

        K = blockdiag(K1, K2)
        M = blockdiag(M1, M2)
    end
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
    I = case.opt.Current(t) / case.param_dim.cell.area / param.scale.I_typ
    j_n = - I / param.NE.as / param.NE.thickness
    j_p = I / param.PE.as / param.PE.thickness
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
    u_n = param.NE.U(cn_surf)
    u_p = param.PE.U(cp_surf)
    j0_n =  param.NE.k * Arrhenius(param.NE.Eac_k, T) * sqrt.(cn_surf .* param.EL.ce0 .* (1.0 .- cn_surf))
    j0_p =  param.PE.k * Arrhenius(param.PE.Eac_k, T) * sqrt.(cp_surf .* param.EL.ce0 .* (1.0 .- cp_surf))
    eta_n = 2 * T * asinh.(j_n / 2.0 ./ j0_n)
    eta_p = 2 * T * asinh.(j_p / 2.0 ./ j0_p)
    V_cell = u_p - u_n + eta_p - eta_n 
    variables["cell voltage"] = V_cell 
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
    return variables
end