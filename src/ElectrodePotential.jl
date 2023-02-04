function  ElectrodePotential(electrode::Electrode, mesh::Mesh, mlen::Int64, T::Float64=1.0)
    """
        Generate the system equatiuon of M du/dt=Ku+F
        inputs = electrode -- eledtrode parameter
                 mesh -- mesh
                 mlen -- length of unknowns
                 T -- temperature
        outputs = M -- mass matrix
                  K -- stiffness matrix
    """
    M = spzeros(mlen, mlen)

    Vi = mesh.element[mesh.gs.ele,:]
    Vj = mesh.element[mesh.gs.ele,:]

    sig_eff = electrode.sig * electrode.eps_s
    coeff = sig_eff .* mesh.gs.weight .* mesh.gs.detJ
    K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff , mlen)
    return M, K
end