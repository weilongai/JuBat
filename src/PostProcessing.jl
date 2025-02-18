function PostProcessing(case::Case, variables::Dict{String, Union{Array{Float64},Float64}}, v::Int64)
    result = Dict()
    result["time [s]"]= variables["time"][1:v] * case.param.scale.t0 
    result["cell voltage [V]"]= variables["cell voltage"][1:v] * case.param.scale.phi
    result["cell current [A]"]= variables["cell current"][1:v] * case.param_dim.cell.I1C
    result["temperature [K]"] = variables["temperature"][1:v] * case.param_dim.scale.T_ref
    result["negative particle center radial stress[Pa]"] = variables["negative particle center radial stress"][:,1:v] * case.param.scale.E_n
    result["positive particle center radial stress[Pa]"] = variables["positive particle center radial stress"][:,1:v] * case.param.scale.E_p
    result["negative particle surface tangential stress[Pa]"] = variables["negative particle surface tangential stress"][:,1:v] * case.param.scale.E_n
    result["positive particle surface tangential stress[Pa]"] = variables["positive particle surface tangential stress"][:,1:v] * case.param.scale.E_p
    result["negative particle surface displacement[m]"] = variables["negative particle surface displacement"][:,1:v] * case.param.scale.r0
    result["positive particle surface displacement[m]"] = variables["positive particle surface displacement"][:,1:v] * case.param.scale.r0
    if case.opt.model == "SPM"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][1:v] * case.param.scale.cp_max
    elseif case.opt.model == "SPMe"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][1:v] * case.param.scale.cp_max
        result["negative electrode exchange current density [A/m^2]"] = variables["negative electrode exchange current density"][1:v] * case.param.scale.j
        result["positive electrode exchange current density [A/m^2]"] = variables["positive electrode exchange current density"][1:v] * case.param.scale.j
        result["negative electrode overpotential [V]"]= variables["negative electrode overpotential"][1:v] * case.param.scale.phi
        result["positive electrode overpotential [V]"]= variables["positive electrode overpotential"][1:v] * case.param.scale.phi
        result["electrolyte lithium concentration [mol/m^3]"] = variables["electrolyte lithium concentration"][:,1:v] * case.param.scale.ce
    elseif case.opt.model == "P2D" || case.opt.model == "sP2D"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative electrode exchange current density [A/m^2]"] = variables["negative electrode exchange current density"][:,1:v] * case.param.scale.j
        result["positive electrode exchange current density [A/m^2]"] = variables["positive electrode exchange current density"][:,1:v] * case.param.scale.j
        result["negative electrode overpotential [V]"]= variables["negative electrode overpotential"][:,1:v] * case.param.scale.phi
        result["positive electrode overpotential [V]"]= variables["positive electrode overpotential"][:,1:v] * case.param.scale.phi
        result["electrolyte lithium concentration [mol/m^3]"] = variables["electrolyte lithium concentration"][:,1:v] * case.param.scale.ce
        result["negative electrode potential [V]"] = variables["negative electrode potential"][:,1:v] * case.param.scale.phi
        result["positive electrode potential [V]"] = variables["positive electrode potential"][:,1:v] * case.param.scale.phi
        result["electrolyte potential in negative electrode [V]"] = variables["electrolyte potential in negative electrode"][:,1:v] * case.param.scale.phi
        result["electrolyte potential in positive electrode [V]"] = variables["electrolyte potential in positive electrode"][:,1:v] * case.param.scale.phi
        result["electrolyte potential [V]"] = variables["electrolyte potential"][:,1:v] * case.param.scale.phi
        result["negative electrode open circuit potential [V]"] = variables["negative electrode open circuit potential"] * case.param.scale.phi
        result["positive electrode open circuit potential [V]"] = variables["positive electrode open circuit potential"] * case.param.scale.phi
        result["negative electrode interfacial current density [A/m^2]"]  = variables["negative electrode interfacial current density"] * case.param.scale.j
        result["positive electrode interfacial current density [A/m^2]"]  = variables["positive electrode interfacial current density"] * case.param.scale.j
    end
    return result
end