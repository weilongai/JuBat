function SetCase(param_dim::Params, opt::Option, y0::Array=[])
    """
    A function to set up the case study
        Input: param_dim::Params, opt::Option
        Output: case 
    """

    param = NormaliseParam(param_dim)
    # |--negative--|--separator--|--positive--|
    if opt.model == "SPM" || opt.model == "SPMe"
        # negative particle
        mesh_np = SetMesh([0, param.NE.rs], opt.Nrn, opt.meshType, opt.gsorder)
        # positive particle
        mesh_pp = SetMesh([0, param.PE.rs], opt.Nrp, opt.meshType, opt.gsorder)
        mesh = Dict("negative particle"=> mesh_np, "positive particle"=>mesh_pp)
        index = Dict(
            "negative particle lithium concentration"=> collect(1:mesh_np.nlen),
            "positive particle lithium concentration"=> mesh_np.nlen .+ collect(1:mesh_pp.nlen),
            "negative particle surface lithium concentration"=> mesh_np.nlen,
            "positive particle surface lithium concentration"=> mesh_np.nlen + mesh_pp.nlen,
            )
        v0 = mesh_np.nlen + mesh_pp.nlen
        if opt.model == "SPMe"
            # electrolyte
            space = [0, param.NE.thickness,  param.NE.thickness + param.SP.thickness, param.PE.thickness + param.SP.thickness + param.NE.thickness]
            mesh_el = SetMesh(space, [opt.Nn, opt.Ns, opt.Np], opt.meshType, opt.gsorder)
            
            mesh_el_ne = PickElement(mesh_el, collect(1:opt.Nn))
            mesh_el_sp = PickElement(mesh_el, collect(opt.Nn + 1:opt.Nn + opt.Ns))
            mesh_el_pe = PickElement(mesh_el, collect(opt.Nn + opt.Ns + 1:opt.Nn + opt.Ns +  opt.Np))
            mesh["electrolyte"] = mesh_el
            mesh["negative electrode"] = mesh_el_ne
            mesh["separator"] = mesh_el_sp
            mesh["positive electrode"] = mesh_el_pe
            index["electrolyte lithium concentration"] = v0 .+ collect(1:mesh_el.nlen)
            index["electrolyte lithium concentration in negative electrode"] = v0 .+ collect(1:mesh_el_ne.nlen)
            index["electrolyte lithium concentration in positive electrode"] = v0 + mesh_el.nlen .- collect(mesh_el_pe.nlen - 1:-1:0)
            index["electrolyte lithium concentration in separator"] = v0 + mesh_el_ne.nlen .+ collect(1: mesh_el_sp.nlen)
            v0 += mesh_el.nlen
        end
    elseif opt.model == "P2D"
        # negative particle
        mesh_np = SetMesh([0, param.NE.rs], opt.Nrn, opt.meshType, opt.gsorder)
        # positive particle
        mesh_pp = SetMesh([0, param.PE.rs], opt.Nrp, opt.meshType, opt.gsorder)
        # electrolyte   
        space = [0, param.NE.thickness,  param.NE.thickness + param.SP.thickness, param.PE.thickness + param.SP.thickness + param.NE.thickness]
        mesh_el = SetMesh(space, [opt.Nn, opt.Ns, opt.Np], opt.meshType, opt.gsorder)
        mesh_el_ne = PickElement(mesh_el, collect(1:opt.Nn))
        mesh_el_sp = PickElement(mesh_el, collect(opt.Nn + 1:opt.Nn + opt.Ns))
        mesh_el_pe = PickElement(mesh_el, collect(opt.Nn + opt.Ns + 1:opt.Nn + opt.Ns +  opt.Np))
        mesh_nps = MultipleMesh(mesh_np, mesh_el_ne.nlen)
        mesh_pps = MultipleMesh(mesh_pp, mesh_el_pe.nlen)

        mesh = Dict(
            "negative particle" => mesh_nps, 
            "positive particle" => mesh_pps,
            "electrolyte" => mesh_el,
            "negative electrode" => mesh_el_ne,
            "separator" => mesh_el_sp,
            "positive electrode" => mesh_el_pe,
            )
        index = Dict(
            "negative particle lithium concentration"=> collect(1:mesh_nps.nlen),
            "positive particle lithium concentration"=> mesh_nps.nlen .+ collect(1:mesh_pps.nlen),
            "negative particle surface lithium concentration"=> collect(mesh_np.nlen:mesh_np.nlen:mesh_nps.nlen),
            "positive particle surface lithium concentration"=> mesh_nps.nlen .+ collect(mesh_pp.nlen:mesh_pp.nlen:mesh_pps.nlen),
            )
        v0 = mesh_nps.nlen + mesh_pps.nlen
        index["electrolyte lithium concentration"] = v0 .+ collect(1:mesh_el.nlen)
        index["electrolyte lithium concentration in negative electrode"] = v0 .+ collect(1:mesh_el_ne.nlen)
        index["electrolyte lithium concentration in positive electrode"] = v0 + mesh_el.nlen .- collect(mesh_el_pe.nlen - 1:-1:0)
        index["electrolyte lithium concentration in separator"] = v0 + mesh_el_ne.nlen - 1 .+ collect(1: mesh_el_sp.nlen)
        v0 += mesh_el.nlen
    end 
    if opt.thermalmodel == "lumped"
        index["temperature"] = [v0 + 1]
        v0 += 1
    end
    if opt.model == "P2D"
        index["negative electrode potential"] = v0 .+ collect(1:mesh_el_ne.nlen)
        index["positive electrode potential"] = v0 .+ mesh_el_ne.nlen .+ collect(1:mesh_el_pe.nlen)
        v0 += mesh_el_ne.nlen + mesh_el_pe.nlen
        index["electrolyte potential"] =  v0 .+ collect(1:mesh_el.nlen)
        index["electrolyte potential in negative electrode"] =  v0 .+ collect(1:mesh_el_ne.nlen)
        index["electrolyte potential in positive electrode"] =  v0 .+ mesh_el.nlen .- collect(mesh_el_pe.nlen - 1:-1:0)
        index["electrolyte potential in separator"] =  v0 + mesh_el_ne.nlen - 1 .+ collect(1: mesh_el_sp.nlen)
    end
    case = Case(param_dim, param, opt, mesh, index)
    return case
end


mutable struct Case
    param_dim::Params   # parameters
    param::Params   # dimensionless parameters
    opt::Option # option for solver
    mesh::Dict{String, Mesh}    # mesh for discretisation
    index::Dict{String, Union{Array{Int64}, Int64}} # the index of unknowns
end