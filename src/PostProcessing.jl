function PostProcessing(case::Case, yt::Array{Float64}, variables::Dict{String, Matrix{Float64}}, v::Integer)
    if case.opt.Model == "SPM"
        result = Dict()
        result["negative particle lithium concentration [mol/m^3]"] = yt[case.index["csn"],1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = yt[case.index["csp"],1:v] * case.param.scale.cp_max
        result["time [s]"]= variables["time"][1:v] * case.param.scale.t0 
        result["cell voltage [V]"]= variables["cell voltage"][1,1:v] * case.param.scale.phi
    else
        error("Error: model $(case.opt.Model)  has not been implement!\n" )
    end
    return result
end