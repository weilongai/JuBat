function Arrhenius(Eac::Float64, T::Union{Float64, Array{Float64}})
"""
    Arrhenius function 
        Input: Eac::Float64 - activation energy (normalised with Eac/(R T0))
                T::Float64 - temperature (normalised with T/T0)
        Output: the results of Arrhenius function
"""
    return exp.(Eac * (1 .- 1 ./ T))
end


function IntV(x::Array{Float64}, mesh::Mesh)
"""
    Integration function over a domain (mesh) 
        Input: x::Array{Float64} - the coefficients to integrate
                mesh::Mesh - mesh structure
        Output: the result of (x) over "mesh"
"""
    return sum(x .* mesh.gs.weight .* mesh.gs.detJ)
end