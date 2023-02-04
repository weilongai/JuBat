function  ElectrolytePotential(param::Params, mesh::Mesh, mlen::Int64, variables::Dict{String, Union{Array{Float64},Float64}})
    """
        Generate the system equatiuon of M du/dt=Ku+F
        inputs = electrode -- eledtrode parameter
                 mesh -- mesh
                 mlen -- length of unknowns
                 variables -- dictionary of variables
        outputs = M -- mass matrix
                  K -- stiffness matrix
    """
    M = spzeros(mlen, mlen)
    Vi = mesh.element[mesh.gs.ele,:]
    Vj = mesh.element[mesh.gs.ele,:]
    ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"] 
    ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"] 
    ce_sp_gs = variables["electrolyte lithium concentration at separator Gauss point"]
    coeff_ne = param.EL.kappa(ce_n_gs) * param.NE.eps ^ param.NE.brugg
    coeff_sp = param.EL.kappa(ce_sp_gs) * param.SP.eps ^ param.SP.brugg
    coeff_pe = param.EL.kappa(ce_p_gs) * param.PE.eps ^ param.PE.brugg
    coeff = - [coeff_ne; coeff_sp; coeff_pe] .* mesh.gs.weight .* mesh.gs.detJ
    K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff, mlen)
    return M, K
end