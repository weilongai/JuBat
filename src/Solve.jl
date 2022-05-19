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
    t = t0 + dt  
    yold = y0 
    M, Kold, Fold = CallModel(case, yold, t) 
    Mt = M - theta * Kold * dt 
    Kt = (1 - theta) * Kold * dt + M 
    Ft = Fold * dt 
    num = round(Integer, (t_end - t0)/dt * 2)
    yt = zeros(size(M,1), num) 
    time = zeros(num,1) 
    v = 1  
    yt[:,1] = y0 
    vt = 2 
    print( "start to solve the problem \n")

    # run the model
    while t <= t_end
        if case.opt.Jacobi == "update"
            Mnew, Knew, Fnew = CallModel(case, yold, t) 
            Mt = M - theta * Knew * dt 
            Kt = (1 - theta) * Kold * dt + M 
            Ft = theta * Fnew * dt + (1 - theta) * Fold * dt 
        end
        ynew = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * yold + Ft) 

        # record the results
        if case.opt.OutputType == "auto" || case.opt.OutputTime == []
            v = v + 1 
            yt[:,v] = ynew 
            time[v,1] = t 
        elseif abs(t - case.opt.OutputTime[vt]) < 1e-7
            v = v + 1 
            yt[:,v] = ynew 
            time[v,1] = t 
            vt = vt + 1 
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
    end
 
    yt = yt[:,1:v] 
    time = time[1:v,1] 
    result = PostProcessing(case, yt, time) 
    print("finish the simulation\n") 
    return result
end

function CallModel(case::Case, y::Array{Float64}, t::Float64)
    if case.opt.Model == "SPM"
        M, K, F = SPM(case, y, t) 
    else
        error( "Error: $(case.opt.Model) model has not been implement!\n ")
    end
    return M, K, F
end

function ModelInitial(case::Case)
    if case.opt.Model == "SPM"
        n1 = case.mesh["negative particle"].nlen
        csn100 = case.param.NE.theta_100
        n2 = case.mesh["positive particle"].nlen
        csp100 = case.param.PE.theta_100
        y0 = [ones(Float64, n1, 1) * csn100;  ones(Float64, n2, 1) * csp100]
    else
        error( "Error: $(case.opt.Model{1}) model has not been implement!\n ")
    end
    return y0
end
