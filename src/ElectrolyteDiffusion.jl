function  ElectrolyteDiffusion(param::Params, mesh::Mesh, mlen::Int64, variables::Dict{String, Union{Array{Float64},Float64}})
    """
        Generate the system equatiuon of M du/dt=Ku+F
        inputs = electrode -- eledtrode parameter
                 mesh -- mesh
                 mlen -- length of unknowns
                 variables -- dictionary of variables
        outputs = M -- mass matrix
                  K -- stiffness matrix
    """
    Vi = mesh.element[mesh.gs.ele,:]
    Vj = mesh.element[mesh.gs.ele,:]

    ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"] 
    ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"] 
    ce_sp_gs = variables["electrolyte lithium concentration at separator Gauss point"]
    coeff_ne = ones(size(ce_n_gs)) * param.NE.eps
    coeff_sp = ones(size(ce_sp_gs)) * param.SP.eps
    coeff_pe = ones(size(ce_p_gs)) * param.PE.eps
    coeff = [coeff_ne; coeff_sp; coeff_pe] .* mesh.gs.weight .* mesh.gs.detJ

    M = Assemble(Vi, Vj, mesh.gs.Ni, mesh.gs.Ni, coeff , mlen)
    T = variables["temperature"]
    De_ne =  param.EL.De(ce_n_gs) * param.NE.eps ^ param.NE.brugg
    De_sp =  param.EL.De(ce_sp_gs) * param.SP.eps ^ param.SP.brugg
    De_pe =  param.EL.De(ce_p_gs) * param.PE.eps ^ param.PE.brugg
    De_eff =  [De_ne; De_sp; De_pe] * Arrhenius(param.EL.Eac_D, T)
    coeff = - De_eff .* mesh.gs.weight .* mesh.gs.detJ
    K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff , mlen)
    return M, K
end