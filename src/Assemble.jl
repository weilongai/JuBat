function Assemble(Vi::Array{Integer}, Vj::Array{Integer}, Ni::Array{Float64}, Nj::Array{Float64}, coeff::Array{Float64}, mlen1::Integer, mlen2::Integer=mlen1)
    # A function to assemble the system matrix
    # Inputs = Vi, Vj, Ni, Nj, coeff, mlen1 and mlen2(=mlen1 if not claimed)
    # To get the integration of 'Ni*Nj*coeff*weight*detJ'
    # to obtian the system matrix M(Vi, Vj) 

    gslen = size(Ni,1)
    gslen1 = size(Ni, 2)
    gslen2 = size(Nj, 2)
    KI = vec(zeros(Integer, gslen * gslen1 * gslen2,1))
    KJ = deepcopy(KI)
    KV = vec(zeros(Float64, gslen * gslen1 * gslen2,1))
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
    
    # ## below is for testing, with speed sacrifice
    # K = zeros(mlen1, mlen2)
    # for i = 1 : gslen
    #     K(Vi(i, 1), Vj(i, 1)) = K(Vi(i, 1), Vj(i, 1)) + Ni(i,1) * Nj(i,1) * coeff(i)
    #     K(Vi(i, 1), Vj(i, 2)) = K(Vi(i, 1), Vj(i, 2)) + Ni(i,1) * Nj(i,2) * coeff(i)
    #     K(Vi(i, 2), Vj(i, 1)) = K(Vi(i, 2), Vj(i, 1)) + Ni(i,2) * Nj(i,1) * coeff(i)
    #     K(Vi(i, 2), Vj(i, 2)) = K(Vi(i, 2), Vj(i, 2)) + Ni(i,2) * Nj(i,2) * coeff(i)
    # end
  