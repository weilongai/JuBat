function SetMesh(domain::Any, num::Any, type::String, gsorder::Int64=4)
"""
    A function to set up mesh
     inputs are 'domain, num, type, gsorder'
     'domain' is an array for the problem domain information
     'num' is a vector for the number of elements
     'type' is a string for the element type
     'gsorder' is int for the order of Gauss quadrature, with default value 4
     e.g. mesh =  SetMesh([0,1,2],[5,5],'L2')
     There should be more sophisticated methods to build a mesh, to be finished
"""

    if type in ["L2", "L3"]
        mesh = Mesh1D(domain, num, type, gsorder)
    else
        error("Error: element type $type has not been implemented!\n")
    end
    return mesh
end

function Mesh1D(domain::Vector{Float64}, num::Any, type::String= "L2", gsorder::Int64= 4)
# to build 1D mesh and its gauss points
    dim = 1
    element_number = round(Int, sum(num))
    if type == "L2"
        ele_node = 2
    elseif type == "L3"
        ele_node = 3
    end
    node = zeros(Float64, element_number * (ele_node-1) + 1, dim)
    v = 0
    for i = 1:length(domain) - 1
        dx = 1/num[i]/(ele_node - 1)
        temp = collect(0:dx:1) .* (domain[i+1] - domain[i]) .+ domain[i]
        len = num[i] * (ele_node - 1)
        if i == 1
            node[v + 1:v + len + 1, 1] = temp
            v = v + len + 1
        else
            node[v + 1:v + len, 1] = temp[2:end,1]
            v = v + len
        end
    end
    nlen = deepcopy(v)
    element = zeros(Int64, element_number, ele_node)
    v = 0 
    ele = 0
    for i in num
        for j = 1:i
            ele = ele + 1
            element[ele, 1:ele_node] = v + 1:v + ele_node
            v = v + ele_node - 1
        end
    end
    gs = GetGS(element, node, gsorder, type)
    mesh = Mesh(type, dim, node, nlen, element, gs)
    return mesh
end

function PickElement(mesh::Mesh, v::Vector{Int64})
    """
        pick elements from a mesh
        Inputs = mesh::Mesh
                v::Vector{Inveger}, the index of elements to be picked up
        outputs = mesh1::Mesh
            mesh1 = mesh(v), including the collection of Gaussian points
        the new mesh will resort the index of node, element and Gauss points 
    """
    type = mesh.type
    dim = mesh.dimension
    element = deepcopy(mesh.element[v,:])
    node_pool = sort(unique(reshape(element,:,1)))
    node_pool_pair = zeros(Int64, maximum(node_pool))
    nlen = length(node_pool)
    node_pool_pair[node_pool] = collect(1:nlen)
    node = mesh.node[node_pool,:]
    element = node_pool_pair[element]
    gsorder = mesh.gs.order
    gs = deepcopy(mesh.gs)
    len = gsorder ^ dim
    v_gs = zeros(Int64, length(v) * len)
    gs.ele = zeros(Int64, length(v) * len)
    for i = 1:len
        v_gs[i:len:length(v) * len] = (v .- 1) .* len .+ i
        gs.ele[i:len:length(v) * len] = 1:length(v)
    end
    gs.x = gs.x[v_gs,:]
    gs.xi = gs.xi[v_gs,:]
    gs.weight = gs.weight[v_gs]
    gs.detJ = gs.detJ[v_gs]
    gs.Ni = gs.Ni[v_gs,:]
    gs.dNidx = gs.dNidx[v_gs,:]
    nlen = length(unique(element))
    mesh_picked = Mesh(type, dim, node, nlen, element, gs)
    return mesh_picked
end

function CombineMesh(meshes::Vector{Mesh})
    """
        combine meshes to a big mesh
        Inputs = meshes::Vector{Mesh}
        outputs = mesh_combined::Mesh
           e.g. meshnew = Combine([mesh1, mesh2, mesh3])
    """
    type = meshes[1].type
    dimension = meshes[1].dimension
    gsorder = meshes[1].gs.order
    for i = 2:length(meshes)
        if type != meshes[i].type
            error("the types of meshes do not match for combination!")
        end
        if dimension != meshes[i].dimension
            error("the dimensions of meshes do not match for combination!")
        end
        if gsorder != meshes[i].gs.order
            error("the orders of Gaussian quadrature do not match for combination!")
        end
    end

    n_mesh = length(meshes)
    n_element = zeros(Int64, n_mesh)
    n_node = zeros(Int64, n_mesh)
    n_gs = zeros(Int64, n_mesh)
    len_node = 0
    len_element = 0
    len_gs = 0
    for i = 1:n_mesh
        n_element[i] = size(meshes[i].element, 1)
        n_node[i] = meshes[i].nlen
        n_gs[i] = size(meshes[i].gs.x, 1)
        len_element += n_element[i]
        len_node += n_node[i]
        len_gs += n_gs[i]
    end

    element = zeros(Int64, len_element, size(meshes[1].element,2))
    node = zeros(Float64, len_node, size(meshes[1].node,2))
    x = zeros(Float64, len_gs, size(meshes[1].gs.x, 2))
    xi = zeros(Float64, len_gs, size(meshes[1].gs.xi, 2))
    weight = zeros(Float64, len_gs)
    detJ = zeros(Float64, len_gs)
    ele = zeros(Int64, len_gs)
    Ni = zeros(Float64, len_gs, size(meshes[1].gs.Ni, 2))
    dNidx = zeros(Float64, len_gs, size(meshes[1].gs.dNidx, 2))

    v_ele = 0
    v_node = 0
    v_gs = 0
    for i = 1:n_mesh
        element[v_ele + 1:v_ele + n_element[i],:] = meshes[i].element .+ v_node
        ele[v_gs + 1:v_gs + n_gs[i]]  = meshes[i].gs.ele .+ v_ele 
        node[v_node + 1:v_node + n_node[i],:] = meshes[i].node
        x[v_gs + 1:v_gs + n_gs[i], :] = meshes[i].gs.x
        xi[v_gs + 1:v_gs + n_gs[i], :]  = meshes[i].gs.xi
        weight[v_gs + 1:v_gs + n_gs[i]]  = meshes[i].gs.weight
        detJ[v_gs + 1:v_gs + n_gs[i]]  = meshes[i].gs.detJ
        Ni[v_gs + 1:v_gs + n_gs[i], :]  = meshes[i].gs.Ni
        dNidx[v_gs + 1:v_gs + n_gs[i], :]  = meshes[i].gs.dNidx
        v_ele += n_element[i]
        v_node += n_node[i]
        v_gs += n_gs[i]
    end

    gs = GaussPoint(x, xi, weight, detJ, ele, Ni, dNidx, gsorder)
    mesh_combined = Mesh(type, dimension, node, len_node, element, gs)
    return mesh_combined
end

function MultipleMesh(mesh::Mesh, n::Int64)
    """
        A function to duplicate a mesh by n times and combine them to one big mesh
            Inputs = mesh::Mesh
                    n::Int64
            Outputs = meshnew::Mesh
    """
    meshes = [mesh]
    for i = 1:n-1
        push!(meshes, mesh)
    end
    meshnew = CombineMesh(meshes)
    return meshnew
end

function GetGS(element::Array{Int64}, node::Array{Float64}, order::Int64, type::String, v=collect(1:size(element,1)))
    if type == "L2" 
        dimen=1
        points = 1:2
    elseif type == "L3"
        dimen=1
        points = [1, 3]
    elseif type == "Q4"
        dimen=2
        points = 1:4
    elseif type == "B8"
        dimen=3
        points = 1:8
    end
    total_num = size(element,1) * order ^ dimen
    x = zeros(Float64, total_num ,dimen)
    weight = zeros(Float64, total_num)
    detJ = zeros(Float64, total_num) 
    ele = zeros(Int64, total_num)
    xi = zeros(Float64, total_num, dimen)
    w, q = GSweight(order,dimen)
    count0 = 0
    for e = 1:size(element, 1)
        sctr = element[e, points]
        for i in eachindex(w)
            pt = q[i, :]
            N, dNdxi = LagrangeBasis(type, dimen, pt)
            J0 = dNdxi * node[sctr, 1:dimen]
            count0 = count0 + 1
            x[count0, 1:dimen] = N * node[sctr, 1:dimen]
            weight[count0] = w[i]
            detJ[count0] = det(J0)
            ele[count0] = v[e]
            xi[count0, 1:dimen] = pt
        end
    end
    Ni, dNi = ShapeFunction1D(element, type, node, xi, ele)
    gs = GaussPoint(x, xi, weight, detJ, ele, Ni, dNi, order)
    return gs
end

function LagrangeBasis(type::String, dimen::Int64, coord::Array{Float64})
    N = zeros(Float64,2^dimen, 1)
    dNdxi = zeros(Float64, 2^dimen, dimen)
    if type == "L2" ||  type == "L3" 
        # 1------2 L2 TWO NODE LINE ELEMENT
        xi = coord[1]
        N[1,1] = (1.0 - xi)/2.0
        N[2,1] = (1.0 + xi)/2.0 
        dNdxi[1,1] = -1.0 /2.0  
        dNdxi[2,1] = 1.0 /2.0 
    elseif type == "Q4"
        ## Q4 FOUR NODE QURARILATERIAL ELEMENT
        # 4---3
        # |   |
        # 1---2
        xi = coord[1] 
        eta = coord[2]
        N=1/4 * [ (1-xi)*(1-eta);
            (1+xi)*(1-eta);
            (1+xi)*(1+eta);
            (1-xi)*(1+eta)]
        dNdxi=1/4 * [-(1-eta)   -(1-xi); 1-eta  -(1+xi); 1+eta  1+xi;   -(1+eta)    1-xi]
    elseif type == "B8"
        ## B4 EIGHT NODE BRICK ELEMENT
        # 4---3   8---7
        # |   |   |   |
        # 1---2 , 5---6
        xi=coord[1]
        eta=coord[2]
        zeta=coord[3]
        I1=1/2 - coord/2 
        I2=1/2 + coord/2
        N=[   I1[1]*I1[2]*I1[3];
            I2[1]*I1[2]*I1[3];
            I2[1]*I2[2]*I1[3];
            I1[1]*I2[2]*I1[3];
            I1[1]*I1[2]*I2[3];
            I2[1]*I1[2]*I2[3];
            I2[1]*I2[2]*I2[3];
            I1[1]*I2[2]*I2[3]   ]
        dNdxi=[-1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
            1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
            1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
            -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;
            -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
            1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
            1+eta+zeta+eta*zeta   1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
            -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ]/8
    end
    
    N=N'
    dNdxi=dNdxi'
    return N, dNdxi
end


function GSweight(order::Int64, dimen::Int64)
    if (order>10 || order<0)
        disp("Order of quadrature too high for Gaussian Quadrature")
    end
    r1pt = zeros(Float64, order) 
    r1wt = zeros(Float64, order)
    W = zeros(Float64, order^dimen)
    Q = zeros(Float64, order^dimen,dimen)
    if order == 1
        r1pt[1] = 0.000000000000000
        r1wt[1] = 2.000000000000000
        
    elseif order == 2
        r1pt[1] = 0.577350269189626
        r1pt[2] =-0.577350269189626
        
        r1wt[1] = 1.000000000000000
        r1wt[2] = 1.000000000000000
        
    elseif order ==  3
        r1pt[1] = 0.774596669241483
        r1pt[2] =-0.774596669241483
        r1pt[3] = 0.000000000000000
        
        r1wt[1] = 0.555555555555556
        r1wt[2] = 0.555555555555556
        r1wt[3] = 0.888888888888889
        
    elseif order ==  4
        r1pt[1] = 0.861134311594053
        r1pt[2] =-0.861134311594053
        r1pt[3] = 0.339981043584856
        r1pt[4] =-0.339981043584856
        
        r1wt[1] = 0.347854845137454
        r1wt[2] = 0.347854845137454
        r1wt[3] = 0.652145154862546
        r1wt[4] = 0.652145154862546
            
    elseif order ==  5
        r1pt[1] = 0.906179845938664
        r1pt[2] =-0.906179845938664
        r1pt[3] = 0.538469310105683
        r1pt[4] =-0.538469310105683
        r1pt[5] = 0.000000000000000
        
        r1wt[1] = 0.236926885056189
        r1wt[2] = 0.236926885056189
        r1wt[3] = 0.478628670499366
        r1wt[4] = 0.478628670499366
        r1wt[5] = 0.568888888888889
            
    elseif order ==  6
        r1pt[1] = 0.932469514203152
        r1pt[2] =-0.932469514203152
        r1pt[3] = 0.661209386466265
        r1pt[4] =-0.661209386466265
        r1pt[5] = 0.238619186003152
        r1pt[6] =-0.238619186003152
        
        r1wt[1] = 0.171324492379170
        r1wt[2] = 0.171324492379170
        r1wt[3] = 0.360761573048139
        r1wt[4] = 0.360761573048139
        r1wt[5] = 0.467913934572691
        r1wt[6] = 0.467913934572691
            
    elseif order ==  7
        r1pt[1] =  0.949107912342759
        r1pt[2] = -0.949107912342759
        r1pt[3] =  0.741531185599394
        r1pt[4] = -0.741531185599394
        r1pt[5] =  0.405845151377397
        r1pt[6] = -0.405845151377397
        r1pt[7] =  0.000000000000000
        
        r1wt[1] = 0.129484966168870
        r1wt[2] = 0.129484966168870
        r1wt[3] = 0.279705391489277
        r1wt[4] = 0.279705391489277
        r1wt[5] = 0.381830050505119
        r1wt[6] = 0.381830050505119
        r1wt[7] = 0.417959183673469
            
    elseif order ==  8
        r1pt[1] =  0.960289856497536
        r1pt[2] = -0.960289856497536
        r1pt[3] =  0.796666477413627
        r1pt[4] = -0.796666477413627
        r1pt[5] =  0.525532409916329
        r1pt[6] = -0.525532409916329
        r1pt[7] =  0.183434642495650
        r1pt[8] = -0.183434642495650
        
        r1wt[1] = 0.101228536290376
        r1wt[2] = 0.101228536290376
        r1wt[3] = 0.222381034453374
        r1wt[4] = 0.222381034453374
        r1wt[5] = 0.313706645877887
        r1wt[6] = 0.313706645877887
        r1wt[7] = 0.362683783378362
        r1wt[8] = 0.362683783378362
        
    else
        r1pt[1] = 0.9739065285
        r1pt[2] = -0.9739065285
        r1pt[3] =  0.8650633677
        r1pt[4] = -0.8650633677
        r1pt[5] =  0.6794095683
        r1pt[6] = -0.6794095683
        r1pt[7] = 0.4333953941
        r1pt[8] =-0.4333953941
        r1pt[9] = 0.1488743390
        r1pt[10] =-0.1488743390
        
        r1wt[1] =0.0666713443
        r1wt[2] =0.0666713443
        r1wt[3] =0.1494513492
        r1wt[4] =0.1494513492
        r1wt[5] =0.2190863625
        r1wt[6] =0.2190863625
        r1wt[7] = 0.2692667193
        r1wt[8] =0.2692667193
        r1wt[9] = 0.2955242247
        r1wt[10] = 0.2955242247
    end
    num=1
    if dimen == 1
        for i = 1:order
            Q[num, 1] = r1pt[i]
            W[num] = r1wt[i]
            num = num+1
        end
    elseif dimen == 2
        for i = 1:order
            for j = 1:order
                Q[num,1] = r1pt[i]   
                Q[num,2] = r1pt[j]
                W[num] = r1wt[i] * r1wt[j]
                num = num + 1
            end
        end
    else
        for i=1:order
            for j=1:order
                for k=1:order
                    Q[num,1] = r1pt[i]
                    Q[num,2] = r1pt[j]
                    Q[num,3] = r1pt[k]
                    W[num] = r1wt[i] * r1wt[j] * r1wt[k]
                    num = num + 1
                end
            end
        end
    end
    return W, Q
end

function ShapeFunction1D(element::Matrix{Int64}, type::String, node::Matrix{Float64}, xi::Array{Float64}, v::Vector{Int64})
    if type == "L3"
            # f1 = x-> (x .- 1).^2 / 4 
            # f2 =  x-> (1 .- x.^2) / 2
            # f3 =  x-> (x .+ 1).^2 / 4
            # df1 = x-> (x .- 1) / 2
            # df2 = x-> -x
            # df3 = x-> (x .+ 1) / 2
            ## another group of shape functions
            f1 = x-> 0.5 * x.^2 - 0.5 * x
            f2 =  x-> - x.^2 .+ 1
            f3 =  x-> 0.5 * x.^2 + 0.5 * x
            df1 = x-> x .- 1 / 2
            df2 = x-> -2 * x
            df3 = x-> x .+ 1 / 2
            
            ele_length = abs.(node[element[v, 3]] - node[element[v, 1]])
            Ni = cat(f1(xi), f2(xi), f3(xi),dims=2)
            dNidX = cat(df1(xi), df2(xi), df3(xi), dims=2)
            dXdx = 2 ./ ele_length * ones(1, 3)
            dNidx = dNidX .* dXdx
    elseif type == "L2"
            f1 =  x-> (1 .- x)/2 
            f2 =  x-> (1 .+ x)/2
            df1 =  x-> -0.5 * ones(Float64,size(x))
            df2 =  x-> 0.5 * ones(Float64,size(x))
            
            ele_length = abs.(node[element[v, 2]] - node[element[v, 1]])
            Ni = cat(f1(xi), f2(xi), dims=2)
            dNidX = cat(df1(xi), df2(xi), dims=2)
            dXdx = 2 ./ ele_length * ones(1, 2)
            dNidx = dNidX .* dXdx
    else
            error("Error: element type $(mesh.type) has not been implemented!\n")
    end
    return Ni, dNidx 
end