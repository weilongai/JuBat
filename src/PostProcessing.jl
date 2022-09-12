function PostProcessing(case::Case, variables::Dict{String, Union{Array{Float64},Float64}}, v::Int64)
    result = Dict()
    result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][1:v] * case.param.scale.cn_max
    result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][1:v] * case.param.scale.cp_max
    result["time [s]"]= variables["time"][1:v] * case.param.scale.t0 
    result["cell voltage [V]"]= variables["cell voltage"][1:v] * case.param.scale.phi

    if case.opt.model == "SPMe"
        result["negative electrode exchange current density [A/m^2]"] = variables["negative electrode exchange current density"][1:v] * case.param.scale.j
        result["positive electrode exchange current density [A/m^2]"] = variables["positive electrode exchange current density"][1:v] * case.param.scale.j
        result["negative electrode overpotential [V]"]= variables["negative electrode overpotential"][1:v] * case.param.scale.phi
        result["positive electrode overpotential [V]"]= variables["positive electrode overpotential"][1:v] * case.param.scale.phi
        result["electrolyte lithium concentration [mol/m^3]"] = variables["electrolyte lithium concentration"][:,v] * case.param.scale.ce
    end
    return result
end