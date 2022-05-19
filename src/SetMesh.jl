mutable struct GaussPoint
    x::Array{Float64}
    xloc::Array{Float64}
    weight::Array{Float64}
    detJ::Array{Float64}
    ele::Array{Integer}
    Ni::Array{Float64}
    dNidx::Array{Float64}
end

mutable struct Mesh
    type::String
    dimension::Integer
    node::Array{Float64}
    nlen::Integer
    element::Array{Integer}
    gs::GaussPoint
end

function SetMesh(domain::Any, num::Any, type::String, gsorder::Integer=4)
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

    if type == "L2"
            ele_node = 2
            mesh = Mesh1D(domain, num, type, gsorder, ele_node)
    elseif type == "L3"
            ele_node = 3
            mesh = Mesh1D(domain, num, type, gsorder, ele_node)
    else
            error("Error: element type $type has not been implemented!\n")
    end
    return mesh
end
    

function Mesh1D(domain, num, type::String= "L2", gsorder::Integer= 4, ele_node::Integer= 2)
# to build 1D mesh and its gauss points
    dim = 1
    element_number = round(Int, sum(num))
    node = zeros(element_number * (ele_node-1) + 1, dim)
    type = type
    dimension = dim
    v = 0
    for i = 1:size(domain,1)
        for j = 1:size(domain,2) - 1
            dx = 1/num[i, j]/(ele_node - 1)
            temp = collect(0:dx:1) .* (domain[1,j+1] - domain[1,j]) .+ domain[1,j]
            len = num[i, j] * (ele_node - 1)
            if j == 1
                node[v + 1:v + len + 1, 1] = temp
                v = v + len + 1
            else
                node[v + 1:v + len, 1] = temp[2:end,1]
                v = v + len
            end
    
        end
    end
    nlen = deepcopy(v)
    element = zeros(Integer, element_number, ele_node)
    v = 0 
    ele = 0
    for i = 1: size(num,1)
        for j = 1:size(num,2)
            for k = 1:num[i,j]
                ele = ele + 1
                element[ele, 1:ele_node] = v + 1:v + ele_node
                v = v + ele_node - 1
            end
        end
        v = v + 1
    end
    gs = GetGS(element[:,[1,ele_node]], node, gsorder, dim)
    mesh = Mesh(type, dim, node, nlen, element, gs)
    return mesh
end


function GetGS(element::Array{Integer}, node::Array{Float64}, order::Integer, dimen::Integer, v=collect(1:size(element,1)))
    total_num = size(element,1) * order ^ dimen
    x = zeros(Float64, total_num ,dimen)
    weight = zeros(Float64, total_num, 1)
    detJ = zeros(Float64, total_num, 1) 
    ele = zeros(Integer, total_num, 1)
    xloc = zeros(Float64, total_num, dimen)
    if dimen==1
        type="L2"
        elen=2
    elseif dimen==2
        type="Q4"
        elen=4
    elseif dimen==3
        type = "B8"
        elen = 8
    end
    w, q = GSweight(order,dimen)
    count0 = 0
    for e = 1:size(element, 1)
        sctr = element[e, 1:elen]
        for i = 1:size(w, 1)
            pt = q[i, :]
            N, dNdxi = LagrangeBasis(type, dimen, pt)
            J0 = dNdxi * node[sctr, 1:dimen]
            count0 = count0 + 1
            x[count0,1:dimen] = N * node[sctr, 1:dimen]
            weight[count0] = w[i]
            detJ[count0] = det(J0)
            ele[count0] = v[e]
            xloc[count0,1:dimen] = pt
        end
    end
    Ni, dNi = ShapeFunction1D(element, type, node, xloc, ele)
    gs = GaussPoint(x, xloc, weight, detJ, ele, Ni, dNi)
    return gs
end

function LagrangeBasis(type::String, dimen::Integer, coord::Array{Float64})
    N = zeros(Float64,2^dimen, 1)
    dNdxi = zeros(Float64, 2^dimen, dimen)
    if type == "L2"
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


function GSweight(order::Integer, dimen::Integer)
    if (order>10 || order<0)
        disp("Order of quadrature too high for Gaussian Quadrature")
    end
    r1pt = zeros(order,1) 
    r1wt = zeros(order,1)
    W = zeros(order^dimen,1)
    Q = zeros(order^dimen,dimen)
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

function ShapeFunction1D(element::Matrix{Integer}, type::String, node::Matrix{Float64}, xloc::Matrix{Float64}, v::Matrix{Integer})
    if type == "L3"
            f1 = x-> (x .- 1).^2 / 4 
            f2 =  x-> (1 .- x.^2) / 2
            f3 =  x-> (x .+ 1).^2 / 4
            df1 = x-> (x .- 1) / 2
            df2 = x-> -x
            df3 = x-> (x .+ 1) / 2
            
            ele_length = node[element[v, 3]] - node[element[v, 1]]
            Ni = cat(f1(xloc), f2(xloc), f3(xloc),dims=2)
            dNidX = cat(df1(xloc), df2(xloc), df3(xloc), dims=2)
            dXdx = 2 ./ ele_length * ones(1, 3)
            dNidx = dNidX .* dXdx
    elseif type == "L2"
            f1 =  x-> (1 .- x)/2 
            f2 =  x-> (1 .+ x)/2
            df1 =  x-> -0.5 * ones(Float64,size(x))
            df2 =  x-> 0.5 * ones(Float64,size(x))
            
            ele_length = node[element[v, 2]] - node[element[v, 1]]
            Ni = cat(f1(xloc), f2(xloc), dims=2)
            dNidX = cat(df1(xloc), df2(xloc), dims=2)
            dXdx = 2 ./ ele_length * ones(1, 2)
            dNidx = dNidX .* dXdx
    else
            error("Error: element type $(mesh.type) has not been implemented!\n")
    end
    return Ni, dNidx 
end