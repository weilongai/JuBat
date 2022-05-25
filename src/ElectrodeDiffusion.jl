function  ElectrodeDiffusion(electrode::Electrode, mesh::Mesh, flux::Array{Float64}, mlen::Integer)
    """
        Generate the system equatiuon of Ma=Ku+F
        where a=du/dt
        inputs = param, mesh, flux, mlen
    """
    Vi = mesh.element[mesh.gs.ele[:,1],:];
    Vj = mesh.element[mesh.gs.ele[:,1],:];
    coeff = mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ;
    M = Assemble(Vi, Vj, mesh.gs.Ni, mesh.gs.Ni, coeff , mlen);
    coeff = - electrode.Ds * mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ;
    K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff, mlen); 
    F =  flux * electrode.Rs ^ 2;
    return [M, K, F]
end