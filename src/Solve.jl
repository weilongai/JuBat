function Solve(case::Case)
    dt_min = case.opt.dt[1] / case.param.scale.t0 
    dt_max = case.opt.dt[2] / case.param.scale.t0 
    RunTime = case.opt.time / case.param.scale.t0 
    t0 = RunTime[1] 
    t_end = RunTime[end]
    
    # initialisation 
    y0 = ModelInitialisation(case) 

    if case.opt.solveType == "Crank-Nicolson"
        theta = 0.5 
    elseif case.opt.solveType == "forward"
        theta = 0 
    elseif "backward"
        theta = 1 
    else
        error( "Error: $(opt.solve_type) difference scheme has not been implemented!\n ") 
    end

    dt = deepcopy(dt_min)
    ddt = deepcopy(dt_min)
    num = round(Int64, (t_end - t0)/dt * 1.5) 
    variables_hist = StandardVariables(case, num)


    t = t0
    vt = 1  
    v = 1 
    yold = deepcopy(y0)
    M, Kold, Fold, variables= CallModel(case, y0, t, jacobi="update") 
    Variable_update!(variables_hist, variables, v) 
    Mt = M - Kold * dt 
    Kt =  M 
    Ft = Fold * dt 
    ynew = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * yold + Ft) 
    t += dt
    print( "start to solve the problem \n")

    # run the model
    while t <= t_end
        if case.opt.jacobi == "update"
            _, Knew, Fnew, variables = CallModel(case, ynew, t, jacobi="update") 
        else
            _, _, Fnew, variables = CallModel(case, ynew, t, jacobi="constant") 
            Knew = deepcopy(Kold)
        end
        Mt = M - theta * Knew * dt 
        Kt = (1 - theta) * Kold * dt + M 
        Ft = theta * Fnew * dt + (1 - theta) * Fold * dt 
        yold = deepcopy(ynew)
        ynew = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * yold + Ft) 

        # record the results
        if case.opt.outputType == "auto" || case.opt.outputTime == [] || abs(t - case.opt.outputTime[vt]) < 1e-7
            v = v + 1 
            Variable_update!(variables_hist, variables, v) 
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
        t += dt 
        if case.opt.jacobi == "update"
            Kold = deepcopy(Knew) 
        end
    end

    result = PostProcessing(case, variables_hist, v) 
    print("finish the simulation\n") 
    return result
end

function CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String="update")
    if case.opt.model == "SPM"
        M, K, F, variables = SPM(case, yt, t, jacobi=jacobi) 
    elseif case.opt.model == "SPMe"
        M, K, F, variables = SPMe(case, yt, t, jacobi=jacobi)
    elseif case.opt.model == "P2D"
        M, K, F, variables = P2D(case, yt, t, jacobi=jacobi)     
    else
        error( "Error: $(case.opt.model) model has not been implemented!\n ")
    end
    return M, K, F,variables
end


