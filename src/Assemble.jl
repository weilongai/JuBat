function Assemble(Vi::Array{Int64}, Vj::Array{Int64}, Ni::Array{Float64}, Nj::Array{Float64}, coeff::Array{Float64}, mlen1::Int64, mlen2::Int64=mlen1)
"""
    A function to assemble the system matrix
    To get the integration of 'K(Vi, Vj) = Ni * Nj * coeff * weight * detJ'
    to obtian the   
    Inputs = Vi, Vj, Ni, Nj, coeff, mlen1 and mlen2(=mlen1 if not claimed)
    Outputs = K -- system matrix K(Vi, Vj)   
"""
    gslen = size(Ni,1)
    gslen1 = size(Ni, 2)
    gslen2 = size(Nj, 2)
    KI = zeros(Int64, gslen * gslen1 * gslen2)
    KJ = deepcopy(KI)
    KV = zeros(Float64, gslen * gslen1 * gslen2)
    v = 0
    for i = 1 : gslen1
        for j = 1 : gslen2
            KI[v+1:v+gslen] = Vi[:, i]
            KJ[v+1:v+gslen] = Vj[:, j]
            KV[v+1:v+gslen] = Ni[:,i] .* Nj[:,j] .* coeff 
            v = v + gslen
        end
    end
    K = sparse(KI, KJ, KV, mlen1, mlen2)
    return K
end

function Assemble1D(Vi::Array{Int64}, Ni::Array{Float64}, coeff::Array{Float64}, mlen1::Int64)
"""
    A function to assemble the loading vector
    To get the integration of 'F(Vi) = Ni * coeff * weight * detJ'        
    Inputs = Vi, Ni, coeff and mlen1
    Outputs = F -- loading vector F(Vi) 
"""
    gslen1 = size(Ni, 1)
    F = zeros(Float64, mlen1)
    for i = 1 : gslen1
            F[Vi[i, :]] += Ni[i,:] * coeff[i]
    end
    return F
end