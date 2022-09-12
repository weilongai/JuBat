"""
    LG M50 cell parameters

    from the paper

    Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W. Dhammika Widanage, and Emma Kendrick. ["Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion param_dim Models."](https://iopscience.iop.org/article/10.1149/1945-7111/ab9050) Journal of the Electrochemical Society 167 (2020): 080534
    
    and references therein.
"""

# Positive Electrode
PE = Electrode()
PE.theta_100 = 0.263849
PE.theta_0 = 0.853974 
PE.thickness = 75.6e-6  
PE.lambda = 2.1
PE.Ds = 4e-15
PE.rho = 3262
PE.heat_Q = 700
PE.eps = 0.335
PE.eps_fi = 0
PE.brugg = 1.5
PE.k = 3.42e-6
PE.cs_max = 63104
PE.cs_0 = 17038
PE.rs = 5.22e-6
PE.sig = 0.18
PE.Eac_D = 0
PE.Eac_k = 0
PE.alpha = 0.5
PE.U = x-> -0.8090*x .+ 4.4875 - 0.0428*tanh.(18.5138*(x .- 0.5542)) - 17.7326*tanh.(15.7890*(x .- 0.3117)) + 17.5842*tanh.(15.9308*(x .- 0.3120))
PE.dUdT = x-> 279.9800214*(tanh.(15.789*x .- 4.9214313)^2) + 0.79239064*(tanh.(18.5138*x .- 10.26034796)^2) .- 1.4510386800000594 - (280.13037336*(tanh.(15.9308*x .- 4.9704096)^2))

# Negative Electrode
NE = Electrode()
NE.theta_100 = 0.910612
NE.theta_0 = 0.0263472 
NE.thickness = 85.2e-6  
NE.lambda = 1.7
NE.Ds = 3.3e-14
NE.rho = 1657
NE.heat_Q = 700
NE.eps = 0.25
NE.eps_fi = 0
NE.brugg = 1.5
NE.k = 6.48e-7
NE.cs_max = 33133
NE.cs_0 = 29866
NE.rs = 5.86e-6
NE.sig = 215.
NE.Eac_D = 0.
NE.Eac_k = 0.
NE.alpha = 0.5
NE.U = x-> 1.97938*exp.(-39.3631*x) .+ 0.2482 - 0.0909*tanh.(29.8538*(x .- 0.1234)) - 0.04478*tanh.(14.9159*(x .- 0.2769)) - 0.0205*tanh.(30.4444*(x .- 0.6103))
NE.dUdT = x-> -0.62411*((sech.(18.5802 .- 30.4444*x))^2 + 4.34813*(sech.(3.68396 .- 29.8538*x))^2 + 1.07022*(sech.(4.13021 .- 14.9159*x))^2) - 211.794*exp.(-39.3631*x)

# Electrolyte
EL = Electrolyte()
EL.De = x-> 1.769e-10
EL.kappa = x-> 0.1297 * x.^3 - 2.51 *x.^1.5 + 3.329 * x
EL.dlnf_dlnc = x-> 1
EL.rho = 1290
EL.heat_Q = 0
EL.tplus = 0.2594
EL.ce0 = 1000

# Separator
SP = Separator()
SP.thickness = 12e-6
SP.lambda = 0.16
SP.rho = 397
SP.heat_Q = 700
SP.eps = 0.47
SP.eps_fi = 0.
SP.brugg = 1.5

# Cell
cell = Cell()    
cell.length = 1.58
cell.width = 6.5e-2
# cell.wrapper = 0.
cell.I1C = 5.
cell.no_layers = 1
cell.capacity = 5
cell.cooling_surface = 5.31e-3
cell.v_h = 4.2
cell.v_l = 2.5
cell.volume = 2.42e-5
cell.rho = 0.
cell.alphaT = 0.
cell.heat_Q = 0.
cell.h = 0.
cell.T0 = 298

# Positive Current Collector
PCC = CurrentCollector()
PCC.thickness = 16e-6
PCC.lambda = 237.
PCC.rho = 2700.
PCC.heat_Q = 897.
PCC.sig =3.6914e7

# Negative Current Collector
NCC = CurrentCollector()
NCC.thickness = 12e-6
NCC.lambda = 401.
NCC.rho = 8960.
NCC.heat_Q = 385.
NCC.sig = 5.8411e7

# Tab
tab = Tab()
# tab.length
# tab.width 
# tab.area

# Binder
binder = Binder()
# binder.rho = 

# Scale
scale = Scale()

# assemble to "param_dim"
param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)