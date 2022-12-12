function StandardVariables(case::Case, num::Int64)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        Nn = 1
        Np = 1
    elseif case.opt.model == "P2D"
        Nn = case.mesh["negative electrode"].nlen
        Np = case.mesh["positive electrode"].nlen
    end

    variables = Dict{String, Union{Array{Float64}, Float64}}(
        "negative particle lithium concentration" => zeros(Float64, Nrn, num),
        "positive particle lithium concentration" => zeros(Float64, Nrp, num),
        "negative particle averaged lithium concentration" => zeros(Float64, Nn, num),
        "positive particle averaged lithium concentration" => zeros(Float64, Np, num),
        "negative particle surface lithium concentration" => zeros(Float64, Nn, num),
        "positive particle surface lithium concentration" => zeros(Float64, Np, num),    
        "negative electrode porosity" => zeros(Float64, Nn, num),
        "positive electrode porosity" => zeros(Float64, Np, num),
        "negative electrode temperature" => zeros(Float64, Nn, num),
        "positive electrode temperature" => zeros(Float64, Np, num),  
        "negative electrode exchange current density" => zeros(Float64, Nn, num),
        "positive electrode exchange current density" => zeros(Float64, Np, num), 
        "negative electrode interfacial current density" => zeros(Float64, Nn, num),
        "positive electrode interfacial current density" => zeros(Float64, Np, num), 
        "negative electrode overpotential" => zeros(Float64, Nn, num),
        "positive electrode overpotential" => zeros(Float64, Np, num), 
        "negative electrode open circuit potential" => zeros(Float64, Nn, num),
        "positive electrode open circuit potential" => zeros(Float64, Np, num),
        "cell voltage" => zeros(Float64, 1, num),
        "time" => zeros(Float64, 1, num),      
        "cell current" => zeros(Float64, 1, num),      
    )
    # additional variables for SPMe and P2D
    if case.opt.model == "SPMe" || case.opt.model == "P2D"
        Ne = case.mesh["electrolyte"].nlen
        Ne_n = case.mesh["negative electrode"].nlen
        Ne_p = case.mesh["positive electrode"].nlen
        Ne_sp = case.mesh["separator"].nlen
        Ne_ngs = case.opt.Nn * case.opt.gsorder
        Ne_pgs = case.opt.Np * case.opt.gsorder
        Ne_spgs = case.opt.Ns * case.opt.gsorder
        variables["electrolyte lithium concentration"] = zeros(Float64, Ne, num)
        variables["electrolyte lithium concentration in negative electrode"] = zeros(Float64, Ne_n, num)
        variables["electrolyte lithium concentration in positive electrode"] = zeros(Float64, Ne_p, num)
        variables["electrolyte lithium concentration in separator"] = zeros(Float64, Ne_sp, num)
        variables["electrolyte lithium concentration at negative electrode Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["electrolyte lithium concentration at positive electrode Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["electrolyte lithium concentration at separator Gauss point"] = zeros(Float64, Ne_spgs, num)
    end
    # extra variables for P2D
    if case.opt.model == "P2D"
        variables["negative electrode potential"] = zeros(Float64, Nn, num)
        variables["positive electrode potential"] = zeros(Float64, Np, num)
        variables["electrolyte potential"] = zeros(Float64, Ne, num)
        variables["electrolyte potential in negative electrode"] = zeros(Float64, Ne_n, num)
        variables["electrolyte potential in positive electrode"] = zeros(Float64, Ne_p, num)
        variables["electrolyte potential in separator"] = zeros(Float64, Ne_sp, num)
        variables["negative electrode interfacial current at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode interfacial current at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode open circuit potential at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode open circuit potential at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode overpotential at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode overpotential at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode exchange current density at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode exchange current density at Gauss point"] = zeros(Float64, Ne_pgs, num)
    end

    if "temperature" in collect(keys(case.index))
        variables["temperature"] = zeros(Float64, length(case.index["temperature"]), num)
    else
        variables["temperature"] = zeros(Float64, 1, num)
    end

    return variables
end

function Variable_update!(variables_hist::Dict{String, Union{Array{Float64},Float64}}, variables::Dict{String, Union{Array{Float64},Float64}}, v::Int64)
    var_list = collect(keys(variables))
    for i in var_list
        variables_hist[i][:,v] = collect(variables[i])
    end
    return variables_hist
end
