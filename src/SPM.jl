function SPM(case::Case, y0::Array{Float64}, t::Float64)
    # negative elecrode particle diffusion 
    mesh1 = case.mesh["negative particle"]
    flux1 = zeros(Float64, mesh1.nlen, 1)
    flux1[end] = case.opt.Current(t) / case.param_dim.cell.Total_surface / case.param.scale.j

    M1, K1, F1 = ElectrodeDiffusion(case.param.NE, mesh1, flux1, mesh1.nlen)
    M1 = M1 .* case.param.scale.ts_n

    # positive elecrode particle diffusion
    mesh2 = case.mesh["positive particle"]
    flux2 = zeros(Float64, mesh2.nlen, 1)
    flux2[end] = - flux1[end]
    M2, K2, F2 = ElectrodeDiffusion(case.param.PE, mesh2, flux2, mesh2.nlen)
    M2 = M2 .* case.param.scale.ts_p 
    M = blockdiag(M1, M2)
    K = blockdiag(K1, K2)
    F = [F1; F2]
    # now need to solve Ma=Kc+f
    return M, K, F
end