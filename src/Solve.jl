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
    dt_temp = 0
    dt_temp_flag = false
    num = round(Int64, (t_end - t0)/dt * 1.5) 
    variables_hist = StandardVariables(case, num)


    t = t0
    vt = 2  
    v = 1 
    y_old = deepcopy(y0)
    M_old, K_old, F_old, variables= CallModel(case, y0, t, jacobi="update") 
    Variable_update!(variables_hist, variables, v)
    t += dt 
    if case.opt.jacobi == "constant"
        RecordMatrix!(case, M_old, K_old)    # record system matrix information 
    end

    print( "start to solve the problem \n")

    # run the model
    while t <= t_end
        M_new, K_new, F_new, variables = CallModel(case, y_old, t, jacobi="update") 
        Mt = M_new - theta * K_new * dt 
        Kt = (1 - theta) * K_old * dt + M_new 
        Ft = theta * F_new * dt + (1 - theta) * F_old * dt 
        y_new = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * y_old + Ft) 
        
        error_y = norm(y_new - y_old) / dt * dt_min / norm(y_old)
        if error_y > 2 * case.opt.dtThreshold
            # reduce dt to dt/2 and recaculate y_new
            dt = min(dt/2, dt_min)
            t -= dt
        else 
            # record the results
            if case.opt.outputType == "auto" || abs(t - RunTime[vt]) < 1e-7
                v = v + 1 
                Variable_update!(variables_hist, variables, v) 
                if abs(t - RunTime[vt]) < 1e-7
                    vt = min(vt + 1, length(RunTime)) 
                end
            end
            
            # adjust time incremental step dt
            if  case.opt.dtType == "auto" && dt_temp_flag == false
                if error_y < 0.5 * case.opt.dtThreshold
                    dt = min(dt * 2, dt_max) 
                elseif error_y >= 1.5 * case.opt.dtThreshold
                    dt = deepcopy(dt_min)
                elseif error_y > case.opt.dtThreshold
                    dt = max(dt / 2, dt_min) 
                end
            elseif dt_temp_flag
                dt = deepcopy(dt_temp)
                dt_temp_flag = false    
            end
            if t + dt > RunTime[vt] && t < RunTime[vt]
                dt_temp = deepcopy(dt)
                dt = abs(RunTime[vt] - t) 
                dt_temp_flag = true
            end

            # update system information
            y_old = deepcopy(y_new)
            K_old = deepcopy(K_new)
            F_old = deepcopy(F_new)
            t += dt 
        end
    end
    result = PostProcessing(case, variables_hist, v) 
    print("finish the simulation\n") 
    return result
end

function CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    if case.opt.model == "SPM"
        M, K, F, variables = SPM(case, yt, t, jacobi=jacobi) 
    elseif case.opt.model == "SPMe"
        M, K, F, variables = SPMe(case, yt, t, jacobi=jacobi)
    elseif case.opt.model == "P2D"
        M, K, F, variables = P2D(case, yt, t, jacobi=jacobi)     
    else
        error( "Error: $(case.opt.model) model has not been implemented!\n ")
    end
    return M, K, F, variables
end

function RecordMatrix!(case::Case, M::SparseArrays.SparseMatrixCSC{Float64, Int64}, K::SparseArrays.SparseMatrixCSC{Float64, Int64})
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        l_np= case.mesh["negative particle"].nlen
        l_pp= case.mesh["positive particle"].nlen
        case.param.NE.M_d = M[1:l_np, 1:l_np]
        case.param.NE.K_d = K[1:l_np, 1:l_np]
        case.param.PE.M_d = M[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]
        case.param.PE.K_d = K[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]  
    elseif case.opt.model == "P2D"
        l_np= case.mesh["negative particle"].nlen
        l_pp= case.mesh["positive particle"].nlen
        l_el = case.mesh["electrolyte"].nlen
        l_ne = case.mesh["negative electrode"].nlen
        l_pe = case.mesh["positive electrode"].nlen
        case.param.NE.M_d = M[1:l_np, 1:l_np]
        case.param.NE.K_d = K[1:l_np, 1:l_np]
        case.param.PE.M_d = M[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]
        case.param.PE.K_d = K[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]
        case.param.NE.M_d = M[1:l_np, 1:l_np]
        case.param.NE.K_d = K[1:l_np, 1:l_np]
        case.param.PE.M_d = M[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]
        case.param.PE.K_d = K[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]  
    else
        error( "Error: $(case.opt.model) model has not been implemented!\n ")
    end
    return case
end

