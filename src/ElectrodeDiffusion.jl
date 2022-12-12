function  ElectrodeDiffusion(electrode::Electrode, mesh::Mesh, mlen::Int64, T::Float64=1.0)
    """
        Generate the system equatiuon of M du/dt=Ku+F
        inputs = electrode -- eledtrode parameter
                 mesh -- mesh
                 mlen -- length of unknowns
                 T -- temperature
        outputs = M -- mass matrix
                  K -- stiffness matrix
    """
    Vi = mesh.element[mesh.gs.ele,:]
    Vj = mesh.element[mesh.gs.ele,:]
    coeff = mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ
    M = Assemble(Vi, Vj, mesh.gs.Ni, mesh.gs.Ni, coeff , mlen)
    Ds_eff =  electrode.Ds * Arrhenius(electrode.Eac_D, T)
    coeff = - Ds_eff * mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ
    K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff, mlen) 
    return M, K
end