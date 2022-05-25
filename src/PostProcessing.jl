function PostProcessing(case::Case, yt::Array{Float64}, variables::Dict{String, Any})
    if case.opt.Model == "SPM"
        result = Dict()
        result["negative particle lithium concentration [mol/m^3]"] = yt[case.index["csn"],:] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = yt[case.index["csp"],:] * case.param.scale.cp_max
        result["time [s]"]= variables["time"] * case.param.scale.t0 
    else
        error("Error: model $(case.opt.Model)  has not been implement!\n" )
    end
    return result
end