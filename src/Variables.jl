function StandardVariables(case::Case, num::Int64)
    n1 = case.mesh["negative particle"].nlen
    n2 = case.mesh["positive particle"].nlen

    variables = Dict(
        "negative particle lithium concentration" => zeros(Float64, n1, num),
        "positive particle lithium concentration" => zeros(Float64, n2, num),
        "negative particle averaged lithium concentration" => zeros(Float64, 1, num),
        "positive particle averaged lithium concentration" => zeros(Float64, 1, num),
        "negative particle surface lithium concentration" => zeros(Float64, 1, num),
        "positive particle surface lithium concentration" => zeros(Float64, 1, num),
        "negative electrode potential" => zeros(Float64, 1, num),
        "positive electrode potential" => zeros(Float64, 1, num),        
        "negative electrode porosity" => zeros(Float64, 1, num),
        "positive electrode porosity" => zeros(Float64, 1, num),
        #"separator porosity" => zeros(Float64, 1, num),
        "negative electrode temperature" => zeros(Float64, 1, num),
        "positive electrode temperature" => zeros(Float64, 1, num),  
        "negative electrode exchange current density" => zeros(Float64, 1, num),
        "positive electrode exchange current density" => zeros(Float64, 1, num), 
        "negative electrode overpotential" => zeros(Float64, 1, num),
        "positive electrode overpotential" => zeros(Float64, 1, num), 
        "cell voltage" => zeros(Float64, 1, num),
        "time" => zeros(Float64, 1, num),          
    )
    if case.opt.model == "SPMe" || case.opt.model == "DFN"
        n3 = case.mesh["electrolyte"].nlen
        variables["electrolyte lithium concentration"] = zeros(Float64, n3, num)
    end
    return variables
end

function Variable_update!(case::Case, variables::Dict{String, Matrix{Float64}}, v::Int64, yt::Matrix{Float64}, t::Float64)
    if case.opt.model == "SPM"
        SPM_update!(case, variables, v, yt, t) 
    elseif case.opt.model == "SPMe"
        SPMe_update!(case, variables, v, yt, t) 
    else
        error( "Error: $(case.opt.model) model has not been implemented!\n ")
    end

end
