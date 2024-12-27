function ThermalLumped(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    t = variables["time"]
    T = variables["temperature"][1]
    MT = param.cell.mass * param.cell.heat_Q * ones(1,1)    
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        I_app = variables["cell current"]
        eta_n = variables["negative electrode overpotential"][1]
        eta_p = variables["positive electrode overpotential"][end]
        csn_surf = variables["negative particle surface lithium concentration"][1]
        csp_surf = variables["positive particle surface lithium concentration"][end]
        Q_ohm = 0
        Q_rxn = abs(I_app * (eta_p - eta_n) ) # reaction heat is always positive
        Q_rev = abs(I_app) * T * (param.PE.dUdT(csp_surf) - param.NE.dUdT(csn_surf)) 
    else
        eta_n_gs = variables["negative electrode overpotential at Gauss point"]
        eta_p_gs = variables["positive electrode overpotential at Gauss point"]
        j_p_gs = variables["positive electrode interfacial current at Gauss point"]
        j_n_gs = variables["negative electrode interfacial current at Gauss point"]
        csn_surf_gs = variables["negative particle surface lithium concentration at Gauss point"]
        csp_surf_gs = variables["positive particle surface lithium concentration at Gauss point"]
        mesh_ne = case.mesh["negative electrode"]
        mesh_pe = case.mesh["positive electrode"]
        mesh_el = case.mesh["electrolyte"]
        Q_rxn = (IntV(j_p_gs * param.PE.as .* eta_p_gs, mesh_pe)  + IntV(j_n_gs .* eta_n_gs * param.NE.as, mesh_ne))
        Q_rev = (IntV(j_p_gs * param.PE.as .* param.PE.dUdT(csp_surf_gs) * T, mesh_pe) + IntV(j_n_gs * param.NE.as .* param.NE.dUdT(csn_surf_gs) * T, mesh_ne))
        phis_n = variables["negative electrode potential"]
        phis_p = variables["positive electrode potential"]
        phie = variables["electrolyte potential"]
        dphis_n_dx_gs = sum(mesh_ne.gs.dNidx .* phis_n[mesh_ne.element[mesh_ne.gs.ele,:]], dims=2)
        dphis_p_dx_gs = sum(mesh_pe.gs.dNidx .* phis_p[mesh_pe.element[mesh_pe.gs.ele,:]], dims=2)
        dphie_dx_gs = sum(mesh_el.gs.dNidx .* phie[mesh_el.element[mesh_el.gs.ele,:]], dims=2)
        sig_n_eff =  param.NE.sig * param.NE.eps_s
        sig_p_eff =  param.PE.sig * param.PE.eps_s
        Q_ohm_n = IntV(sig_n_eff * dphis_n_dx_gs .^ 2, mesh_ne) 
        Q_ohm_p = IntV(sig_p_eff * dphis_p_dx_gs .^ 2, mesh_pe)
        ce = variables["electrolyte lithium concentration"]
        ce_gs = sum(mesh_el.gs.Ni .* ce[mesh_el.element[mesh_el.gs.ele,:]], dims=2)
        dcedx_gs = sum(mesh_el.gs.dNidx .* ce[mesh_el.element[mesh_el.gs.ele,:]], dims=2)
        tau_ne = param.NE.eps ^ param.NE.brugg * ones(size(j_n_gs))
        tau_pe = param.PE.eps ^ param.PE.brugg * ones(size(j_p_gs))
        tau_sp = param.SP.eps ^ param.SP.brugg * ones(case.opt.Ns * mesh_el.gs.order, 1)
        tau_el = [tau_ne; tau_sp; tau_pe]
        kappa_e_eff = param.EL.kappa(ce_gs,T) .* tau_el
        Q_ohm_e = IntV(kappa_e_eff .* dphie_dx_gs .^ 2 - 2 * kappa_e_eff * T *  (1 - param.EL.tplus) .* param.EL.dlnf_dlnc(ce_gs) ./ ce_gs .* dcedx_gs .* dphie_dx_gs, mesh_el)
        Q_ohm = Q_ohm_n + Q_ohm_p + Q_ohm_e
    end
    q = param.cell.h * param.cell.cooling_surface .* (T -  param.cell.T_amb)
    return MT, Q_rxn + Q_ohm + Q_rev - q
end