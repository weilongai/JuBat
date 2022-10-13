function ThermalLumped(case::Case, variables::Dict{String, Any}, v::Int64, T::Array{Float64})
    param = case.param
    param_dim = case.param_dim
    I_app = case.opt.Current(t) / case.param_dim.cell.area / param.scale.I_typ
    MT = param.cell.mass * param.cell.heat_Q
    eta_n = variables["negative electrode overpotential"][1,v]
    eta_p = variables["positive electrode overpotential"][1,v]
    csn_surf = variables["negative particle surface lithium concentration"][1,v]
    csp_surf = variables["positive particle surface lithium concentration"][1,v]
    # Q_ohm = 0
    Q_rxn = I_app * (eta_p - eta_n) * param.scale.A_cell
    Q_rev = I_app * T * (param.PE.dUdT(csp_surf) - param.NE.dUdT(csn_surf)) * param.scale.A_cell
    q = param.cell.h * param.cell.cooling_surface * (T -  param.cell.T_amb)
    return MT, Q_rxn + Q_rev - q
end