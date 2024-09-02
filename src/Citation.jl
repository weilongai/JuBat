function Citation(model=[])
    # add reference here
    cits = Dict{String, String}(
        "ai2024a" => "W. Ai, Y. Liu, JuBat: A Julia-based framework for battery modelling using finite element method, SoftwareX. 27 (2024) 101760. https://doi.org/10.1016/j.softx.2024.101760",
        "ai2023" => "W. Ai, Y. Liu, Improving the convergence rate of Newman’s battery model using 2nd order finite element method, J. Energy Storage. 67 (2023) 107512. https://doi.org/10.1016/j.est.2023.107512",
        "ai2024b" => "W. Ai, Y. Liu, sP2D: Simplified pseudo 2D battery model by piecewise sinusoidal/quadratic functions of potential curves, J. Energy Storage. 86 (2024) 111386. https://doi.org/10.1016/j.est.2024.111386",
    ) 
    v = ["ai2024a"]
    for i in model
        # add citation here
        if i == "sP2D"
            v = vcat(v, "ai2024b")
        end
        if i == "L3"
            v = vcat(v, "ai2023")
        end
    end

    v = unique(v)
    print( "Reference \n")
    for i in eachindex(v)
        print("[$i] $(cits[v[i]]) \n")
    end
end