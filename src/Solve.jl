function Solve(case::Case)
    dt_min = case.opt.dt[1] / case.param. scale.t0 
    dt_max = case.opt.dt[2] / case.param. scale.t0 
    RunTime = case.opt.Time / case.param. scale.t0 
    t0 = RunTime[1] 
    t_end = RunTime[end] 
    if isempty(case.opt.y0)
        y0 = ModelInitial(case) 
    else
        y0 = case.opt.y0 
    end
    if case.opt.SolveType == "Crank-Nicolson"
        theta = 0.5 
    elseif case.opt.SolveType == "forward"
        theta = 0 
    elseif "backward"
        theta = 1 
    else
        error( "Error: $(opt.Solve_type) difference scheme has not been implement!\n ") 
    end
    # initialisation
    dt = deepcopy(dt_min)
    ddt = deepcopy(dt_min)
   
    num = round(Integer, (t_end - t0)/dt * 2) 
    v = 1  
    variables = StandardVariables(case, num)
    yold = y0 
    SPM_update!(case, variables, v, y0, t0) 
    M = CallModel(case, "M") 
    Kold= CallModel(case, "K") 
    Fold = CallModel_BC(case, t0) 
  
    yt = zeros(size(M,1), num) 
    yt[:,1] = y0 
    time = zeros(1,num) 
    t = t0 + dt
    vt = 2 
    print( "start to solve the problem \n")

    # run the model
    while t <= t_end
        if case.opt.Jacobi == "update"
            Knew = CallModel(case, "F") 
        else
            Knew = Kold
        end
        Fnew = CallModel_BC(case, t)
        Mt = M - theta * Knew * dt 
        Kt = (1 - theta) * Kold * dt + M 
        Ft = theta * Fnew * dt + (1 - theta) * Fold * dt 

        ynew = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * yold + Ft) 

        # record the results
        if case.opt.OutputType == "auto" || case.opt.OutputTime == [] || abs(t - case.opt.OutputTime[vt]) < 1e-7
            v = v + 1 
            yt[:,v] = ynew 
            time[1, v] = t 
            SPM_update!(case, variables, v, ynew, t) 
            if length(case.opt.OutputTime) > 0 && abs(t - case.opt.OutputTime[vt]) < 1e-7
                vt = vt + 1 
            end
        end
        
        # adjust time incremental step dt
        if case.opt.OutputTime != [] && t + dt > case.opt.OutputTime[vt] && t < case.opt.OutputTime[vt]
            dt = abs(case.opt.OutputTime[vt] - t) 
            vt = min(vt + 1, size(case.opt.OutputTime, 2)) 
		elseif  case.opt.dtType == "auto"
            test = norm(ynew - yold) / norm(yold) 
            if test < case.opt.dtThreshold
		        ddt = ddt*2
                dt = min(ddt, dt_max) 
            elseif test > 4 * case.opt.dtThreshold
                dt = deepcopy(dt_min)
		        ddt = deepcopy(dt_min)
            end
        end

        yold = ynew 
        t = t + dt 
        if case.opt.Jacobi == "update"
            Kold = Knew 
        end
    end
 
    yt = yt[:,1:v] 
    result = PostProcessing(case, yt, variables, v) 
    print("finish the simulation\n") 
    return result
end

function CallModel(case::Case, opt::String)
    if case.opt.Model == "SPM"
        K = SPM(case, opt) 
    else
        error( "Error: $(case.opt.Model) model has not been implement!\n ")
    end
    return K
end

function CallModel_BC(case::Case, t::Float64)
    if case.opt.Model == "SPM"
        F = SPM_BC(case, t) 
    else
        error( "Error: $(case.opt.Model) model has not been implement!\n ")
    end
    return F
end


function CallModel_update!(case::Case, variables::Dict{String, Any}, v::Integer, yt::Array{Float64}, t::Float64)
    if case.opt.Model == "SPM"
        SPM_update!(case, variables, v, yt, t) 
    else
        error( "Error: $(case.opt.Model) model has not been implement!\n ")
    end

end

function ModelInitial(case::Case)
    if case.opt.Model == "SPM"
        n1 = case.mesh["negative particle"].nlen
        csn0 = case.param.NE.cs_0
        n2 = case.mesh["positive particle"].nlen
        csp0 = case.param.PE.cs_0
        y0 = [ones(Float64, n1, 1) * csn0;  ones(Float64, n2, 1) * csp0]
    else
        error( "Error: $(case.opt.Model{1}) model has not been implement!\n ")
    end
    return y0
end


function StandardVariables(case::Case, num::Integer)
    n1 = case.mesh["negative particle"].nlen
    n2 = case.mesh["positive particle"].nlen
    if case.opt.Model == "SPM"
        variables = Dict(
            "negative particle lithium concentration" => zeros(n1, num),
            "positive particle lithium concentration" => zeros(n2, num),
            "negative electrode potential" => zeros(1, num),
            "positive electrode potential" => zeros(1, num),
            "negative particle averaged lithium concentration" => zeros(1, num),
            "positive particle averaged lithium concentration" => zeros(1, num),
            "negative particle surface lithium concentration" => zeros(1, num),
            "positive particle surface lithium concentration" => zeros(1, num),
            "negative electrode porosity" => zeros(1, num),
            "positive electrode porosity" => zeros(1, num),
            #"separator porosity" => zeros(1, num),
            "negative electrode temperature" => zeros(1, num),
            "positive electrode temperature" => zeros(1, num),  
            "negative electrode exchange current density" => zeros(1, num),
            "positive electrode exchange current density" => zeros(1, num), 
            "cell voltage" => zeros(1, num),
            "time" => zeros(1, num),          
        )
    else
        error( "Error: $(case.opt.Model{1}) model has not been implement!\n ")
    end
    return variables
end
