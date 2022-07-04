function  ElectrolyteDiffusion(param::Params, mesh::Mesh, mlen::Int64, opt::Option, T::Float64=1.0)
    """
        Generate the system equatiuon of Ma=Ku+F
        where a=du/dt
        inputs = param, mesh, flux, mlen
        outputs = M or K
    """
    Vi = mesh.element[mesh.gs.ele,:]
    Vj = mesh.element[mesh.gs.ele,:]
    if opt.jacobi_M == "constant"
        M = spzeros(mlen, mlen)
    else
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
    end
    ce = param.EL.ce0 # need to modify later
    De_eff =  param.EL.De(ce) * Arrhenius(param.EL.Eac_D, T)
    coeff = - De_eff .* mesh.gs.weight .* mesh.gs.detJ
    K = Assemble(Vi, Vj, mesh.gs.dNidx, mesh.gs.dNidx, coeff , mlen)
    return M, K
end