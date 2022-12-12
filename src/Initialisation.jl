function ModelInitialisation(case::Case)
    if isempty(case.opt.y0)
        if case.opt.model == "SPM"
            Nrn = case.mesh["negative particle"].nlen
            csn0 = ones(Float64, Nrn, 1) * case.param.NE.cs0
            Nrp = case.mesh["positive particle"].nlen
            csp0 = ones(Float64, Nrp, 1) * case.param.PE.cs0
            y0 = [csn0;  csp0]
        elseif case.opt.model == "SPMe"
            Nrn = case.mesh["negative particle"].nlen
            csn0 = ones(Float64, Nrn, 1) *  case.param.NE.cs0
            Nrp = case.mesh["positive particle"].nlen
            csp0 = ones(Float64, Nrp, 1) *  case.param.PE.cs0
            Ne = case.mesh["electrolyte"].nlen
            ce0 = ones(Float64, Ne, 1) *  case.param.EL.ce0
            y0 = [csn0;  csp0; ce0]
        elseif case.opt.model == "P2D"
            Nrn = case.mesh["negative particle"].nlen
            Nrp = case.mesh["positive particle"].nlen
            Ne = case.mesh["electrolyte"].nlen
            Nn = case.mesh["negative electrode"].nlen
            Np = case.mesh["positive electrode"].nlen
            csn0 = ones(Float64, Nrn, 1) * case.param.NE.cs0
            csp0 = ones(Float64, Nrp, 1) * case.param.PE.cs0
            ce0 = ones(Float64, Ne, 1) * case.param.EL.ce0
            phie0 = - ones(Float64, Ne, 1) * case.param.NE.U(case.param.NE.cs0)
            phis_p =  ones(Float64, Nn, 1) * case.param.PE.U(case.param.PE.cs0) .+ phie0[1] # guessed values are not used
            phis_n = zeros(Float64, Np, 1)

            y0 = [csn0;  csp0; ce0; phis_n; phis_p; phie0]
        else
            error( "Error: $(case.opt.model{1}) model has not been implemented!\n ")
        end
    else
        y0 = case.opt.y0 
    end
    return y0
end
