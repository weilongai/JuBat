function SetCase(param_dim::Params, opt::Option)
    """
    A function to set up the case study
        Input: param_dim::Params, opt::Option
        Output: case 
    """

    param = NormaliseParam(param_dim)
    if opt.model == "SPM"
        # negative particle
        mesh1 = SetMesh([0, param.NE.Rs], opt.Nrn, opt.meshType, opt.gsOrder);
        # positive particle
        mesh2 = SetMesh([0, param.PE.Rs], opt.Nrp, opt.meshType, opt.gsOrder);
        mesh = Dict("negative particle"=> mesh1, "positive particle"=>mesh2)
        index = Dict("csn"=> collect(1:opt.Nrn+1),"csp"=> collect(opt.Nrn + 2:opt.Nrn + opt.Nrp + 2))

    elseif opt.model == "SPMe" || opt.model == "DFN"
        # |--negative--|--separator--|--positive--|
        # negative particle
        mesh1 = SetMesh([0, param.NE.Rs], opt.Nrn, opt.meshType, opt.gsOrder);
        # positive particle
        mesh2 = SetMesh([0, param.PE.Rs], opt.Nrp, opt.meshType, opt.gsOrder);
        # electrolyte
        space = [0, param.NE.thickness,  param.NE.thickness + param.SP.thickness, param.PE.thickness + param.SP.thickness + param.NE.thickness]
        mesh3 = SetMesh(space, [opt.Nn, opt.Ns, opt.Np], opt.meshType, opt.gsOrder)
        
        mesh4 = PickElement(mesh3, collect(1:opt.Nn))
        mesh5 = PickElement(mesh3, collect(opt.Nn + 1:opt.Nn + opt.Ns))
        mesh6 = PickElement(mesh3, collect(opt.Nn + opt.Ns + 1:opt.Nn + opt.Ns +  opt.Np))
        mesh = Dict(
            "negative particle"=>mesh1, 
            "positive particle"=>mesh2,
            "electrolyte"=>mesh3, 
            "negative electrode"=>mesh4,
            "separator"=>mesh5,
            "positive electrode"=>mesh6,
            )
        index = Dict(
            "csn"=> collect(1:mesh1.nlen),
            "csp"=> mesh1.nlen .+ collect(1:mesh2.nlen) ,
            "ce"=>  mesh1.nlen + mesh2.nlen .+ collect(1:mesh3.nlen),
            "cen"=> mesh1.nlen + mesh2.nlen .+ collect(1:mesh4.nlen),
            "cep"=> mesh3.nlen + mesh1.nlen + mesh2.nlen + 1 .- collect(1:mesh6.nlen),
            )        
    end 
    case = Case(param_dim, param, opt, mesh, index)
    return case
end


mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Array{Int64}}
end