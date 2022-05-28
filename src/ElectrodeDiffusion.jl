function  ElectrodeDiffusion(electrode::Electrode, mesh::Mesh, mlen::Integer, opt::String)
    """
        Generate the system equatiuon of Ma=Ku+F
        where a=du/dt
        inputs = param, mesh, flux, mlen
    """
    Vi = mesh.element[mesh.gs.ele[:,1],:];
    Vj = mesh.element[mesh.gs.ele[:,1],:];
    if opt == "M"
        coeff = mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ;
        K = Assemble(Vi, Vj, mesh.gs.Ni, mesh.gs.Ni, coeff , mlen);
    elseif opt == "K"
        coeff = - electrode.Ds * mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ;
        K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff, mlen); 
    end
    return K
end