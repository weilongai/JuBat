function Solve(case::Case)
    dt_min = case.opt.dt[1] / case.param. scale.t0 
    dt_max = case.opt.dt[2] / case.param. scale.t0 
    RunTime = case.opt.time / case.param. scale.t0 
    t0 = RunTime[1] 
    t_end = RunTime[end] 
    if isempty(case.opt.y0)
        y0 = ModelInitial(case) 
    else
        y0 = case.opt.y0 
    end
    if case.opt.solveType == "Crank-Nicolson"
        theta = 0.5 
    elseif case.opt.solveType == "forward"
        theta = 0 
    elseif "backward"
        theta = 1 
    else
        error( "Error: $(opt.solve_type) difference scheme has not been implemented!\n ") 
    end
    # initialisation
    dt = deepcopy(dt_min)
    ddt = deepcopy(dt_min)
   
    num = round(Int64, (t_end - t0)/dt * 2) 
    v = 1  
    variables = StandardVariables(case, num)
    yold = y0 
    Variable_update!(case, variables, v, y0, t0)   
    M = CallModel(case, "M") 
    Kold= CallModel(case, "K") 
    Fold = CallModel_BC(case, t0) 
  
    yt = zeros(Float64, size(M,1), num) 
    yt[:,1] = y0 
    time = zeros(Float64, 1,num) 
    t = t0 + dt
    vt = 2 
    print( "start to solve the problem \n")

    # run the model
    while t <= t_end
        if case.opt.jacobi == "update"
            Knew = CallModel(case, "K") 
        else
            Knew = deepcopy(Kold)
        end
        Fnew = CallModel_BC(case, t)
        Mt = M - theta * Knew * dt 
        Kt = (1 - theta) * Kold * dt + M 
        Ft = theta * Fnew * dt + (1 - theta) * Fold * dt 

        ynew = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * yold + Ft) 

        # record the results
        if case.opt.outputType == "auto" || case.opt.outputTime == [] || abs(t - case.opt.outputTime[vt]) < 1e-7
            v = v + 1 
            yt[:,v] = ynew 
            time[1, v] = t 
            Variable_update!(case, variables, v, ynew, t) 
            if length(case.opt.outputTime) > 0 && abs(t - case.opt.outputTime[vt]) < 1e-7
                vt = vt + 1 
            end
        end
        
        # adjust time incremental step dt
        if case.opt.outputTime != [] && t + dt > case.opt.outputTime[vt] && t < case.opt.outputTime[vt]
            dt = abs(case.opt.outputTime[vt] - t) 
            vt = min(vt + 1, size(case.opt.outputTime, 2)) 
		elseif  case.opt.dtType == "auto"
            change = norm(ynew - yold) / norm(yold) 
            if change < case.opt.dtThreshold
		        ddt = ddt*2
                dt = min(ddt, dt_max) 
            elseif change > 4 * case.opt.dtThreshold
                dt = deepcopy(dt_min)
		        ddt = deepcopy(dt_min)
            else
                ddt = ddt/2
                dt = max(ddt, dt_min)  
            end
        end

        yold = ynew 
        t = t + dt 
        if case.opt.jacobi == "update"
            Kold = deepcopy(Knew) 
        end
    end
 
    yt = yt[:,1:v] 
    result = PostProcessing(case, yt, variables, v) 
    print("finish the simulation\n") 
    return result
end

function CallModel(case::Case, opt::String)
    if case.opt.model == "SPM"
        K = SPM(case, opt) 
    elseif case.opt.model == "SPMe"
        K = SPMe(case, opt)    
    else
        error( "Error: $(case.opt.model) model has not been implemented!\n ")
    end
    return K
end

function CallModel_BC(case::Case, t::Float64)
    if case.opt.model == "SPM"
        F = SPM_BC(case, t) 
    elseif case.opt.model == "SPMe"
        F = SPMe_BC(case, t)
    else
        error( "Error: $(case.opt.model) model has not been implemented!\n ")
    end
    return F
end



function ModelInitial(case::Case)
    if case.opt.model == "SPM"
        n1 = case.mesh["negative particle"].nlen
        csn0 = case.param.NE.cs_0
        n2 = case.mesh["positive particle"].nlen
        csp0 = case.param.PE.cs_0
        y0 = [ones(Float64, n1, 1) * csn0;  ones(Float64, n2, 1) * csp0]
    elseif case.opt.model == "SPMe"
        n1 = case.mesh["negative particle"].nlen
        csn0 = case.param.NE.cs_0
        n2 = case.mesh["positive particle"].nlen
        csp0 = case.param.PE.cs_0
        n3 = case.mesh["electrolyte"].nlen
        ce0 = case.param.EL.ce0
        y0 = [ones(Float64, n1, 1) * csn0;  ones(Float64, n2, 1) * csp0; ones(Float64, n3, 1) * ce0]
    else
        error( "Error: $(case.opt.model{1}) model has not been implemented!\n ")
    end
    return y0
end


