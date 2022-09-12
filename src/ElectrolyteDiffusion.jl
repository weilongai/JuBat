function  ElectrolyteDiffusion(param::Params, mesh::Mesh, mlen::Int64, variables::Dict{String, Union{Array{Float64},Float64}})
    """
        Generate the system equatiuon of Ma=Ku+F
        where a=du/dt
        inputs = param, mesh, flux, mlen
        outputs = M or K
    """
    Vi = mesh.element[mesh.gs.ele,:]
    Vj = mesh.element[mesh.gs.ele,:]
    coeff = mesh.gs.weight .* mesh.gs.detJ
    Nn=10
    Ns=10
    Np=10 # need change later when including variables
    v_ne = collect(1:Nn * mesh.gs.order)
    v_sp = Nn * mesh.gs.order .+ collect(1:Ns)
    v_pe = (Nn + Ns) * mesh.gs.order .+ collect(1:Np)
    coeff[v_ne] .*= param.NE.eps
    coeff[v_sp] .*= param.SP.eps
    coeff[v_pe] .*= param.PE.eps
    M = Assemble(Vi, Vj, mesh.gs.Ni, mesh.gs.Ni, coeff , mlen)

    ce = variables["electrolyte lithium concentration"]
    T = variables["temperature"]
    De_eff =  param.EL.De(ce) * Arrhenius(param.EL.Eac_D, T)
    coeff = - De_eff .* mesh.gs.weight .* mesh.gs.detJ
    K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff , mlen)
    return M, K
end