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
        mlen1 = opt.Nrn
        csn = zeros(Float64, mlen1, 1000);

        # positive particle
        mesh2 = SetMesh([0 param.PE.Rs],opt.Nrp, opt.MeshType, opt.gsOrder);
        mlen2 = opt.Nrp
        csp = zeros(Float64, mlen1, 1000);
        mesh = Dict("negative particle"=> mesh1, "positive particle"=>mesh2)
        vars = Dict("negative particle Li concentration"=> csn, "positive particle Li concentration"=>csp)
        index = Dict("csn"=> collect(1:opt.Nrn),"csp"=> collect(opt.Nrn + 1:opt.Nrn + opt.Nrp))
        case = Case(param_dim, param, opt, mesh, vars, index)
        return case
    end
end


mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    vars::Dict{String, Array{Float64}}
    index::Dict{String, Array{Integer}}
end