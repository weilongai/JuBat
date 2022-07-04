function  ElectrodeDiffusion(electrode::Electrode, mesh::Mesh, mlen::Int64, opt::Option, T::Float64=1.0)
    """
        Generate the system equatiuon of Ma=Ku+F
        where a=du/dt
        inputs = param, mesh, flux, mlen
    """
    Vi = mesh.element[mesh.gs.ele,:]
    Vj = mesh.element[mesh.gs.ele,:]
    if opt.jacobi_M == "constant"
        M = spzeros(mlen, mlen)
    else
        coeff = mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ
        M = Assemble(Vi, Vj, mesh.gs.Ni, mesh.gs.Ni, coeff , mlen)
    end
    Ds_eff =  electrode.Ds * Arrhenius(electrode.Eac_D, T)
    coeff = - Ds_eff * mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ
    K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff, mlen) 
    return M, K
end