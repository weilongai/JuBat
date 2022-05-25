function SetCase(param_dim::Params, opt::Option)
    """
    A function to set up the case study
        Input: param_dim::Params, opt::Option
        Output: case 
    """
    if opt.Model == "SPM"

        param = NormaliseParam(param_dim)
        # negative particle
        mesh1 = SetMesh([0 param.NE.Rs], opt.Nrn, opt.MeshType, opt.gsOrder);

        # positive particle
        mesh2 = SetMesh([0 param.PE.Rs],opt.Nrp, opt.MeshType, opt.gsOrder);
        mesh = Dict("negative particle"=> mesh1, "positive particle"=>mesh2)
        index = Dict("csn"=> collect(1:opt.Nrn+1),"csp"=> collect(opt.Nrn + 2:opt.Nrn + opt.Nrp + 2))
        case = Case(param_dim, param, opt, mesh, index)
        return case
    end
end


mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Array{Integer}}
end