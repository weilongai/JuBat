function Citation(cite=[])
    # add reference here
    papers = Dict{String, String}(
        "ai2024a" => "W. Ai, Y. Liu, JuBat: A Julia-based framework for battery modelling using finite element method, SoftwareX. 27 (2024) 101760. https://doi.org/10.1016/j.softx.2024.101760",
        "ai2023" => "W. Ai, Y. Liu, Improving the convergence rate of Newman’s battery model using 2nd order finite element method, J. Energy Storage. 67 (2023) 107512. https://doi.org/10.1016/j.est.2023.107512",
        "ai2024b" => "W. Ai, Y. Liu, sP2D: Simplified pseudo 2D battery model by piecewise sinusoidal/quadratic functions of potential curves, J. Energy Storage. 86 (2024) 111386. https://doi.org/10.1016/j.est.2024.111386",
    ) 
    cite = unique(cite)
    cite = vcat("ai2024a", cite)
    print( "Reference \n")
    for i in eachindex(cite)
        print("[$i] $(papers[cite[i]]) \n")
    end
end