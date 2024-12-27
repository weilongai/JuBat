"""
    Northrop (LiCoO2-LiC6) cell parameters

    from the paper

    Marcello Torchio, Lalo Magni, R. Bhushan Gopaluni, Richard D. Braatz and Davide M. Raimondo. 
    ["LIONSIMBA: A Matlab Framework Based on a Finite Volume Model Suitable for Li-Ion Battery Design, Simulation, and Control."]
    (https://iopscience.iop.org/article/10.1149/2.0291607jes) 
    Journal of The Electrochemical Society, 163 (7) A1192-A1205 (2016)
    
    and references therein.
"""
# Positive Electrode
PE = Electrode()
PE.theta_100 = 0.49550
PE.theta_0 = 0.99174 
PE.thickness = 80e-6  
PE.lambda = 2.1
PE.Ds = 1e-14
PE.rho = 2500
PE.heat_Q = 700
PE.eps = 0.385
PE.eps_fi = 0.025
PE.brugg = 4
PE.k = 2.334e-11 * 96485.33289 # convert from m^2.5/mol^0.5/s (LIONSIMBA) to Am^2.5/mol^1.5 (present)
PE.cs_max = 51554
PE.cs0 = 29252
PE.rs = 2e-6
PE.sig = 100
PE.Eac_D = 5000
PE.Eac_k = 5000
PE.alpha = 0.5
PE.U = x-> (-4.656 .+ 88.669*x.^2 - 401.119*x.^4 + 342.909*x.^6 - 462.471*x.^8 + 433.434*x.^10) ./ (-1 .+ 18.933*x.^2 - 79.532*x.^4 + 37.311*x.^6 - 73.083*x.^8 + 95.96*x.^10)
PE.dUdT = x->  -0.001 * (0.199521039 .- 0.928373822*x + 1.364550689000003*x.^2 - 0.6115448939999998*x.^3) ./ (1 .- 5.661479886999997*x + 11.47636191*x.^2 - 9.82431213599998*x.^3 + 3.048755063*x.^4)

# Negative Electrode
NE = Electrode()
NE.theta_100 = 0.85510
NE.theta_0 = 0.01429 
NE.thickness = 88e-6  
NE.lambda = 1.7
NE.Ds = 3.9e-14
NE.rho = 2500
NE.heat_Q = 700
NE.eps = 0.485
NE.eps_fi = 0.0326
NE.brugg = 4
NE.k = 5.031e-11 * 96485.33289 # convert from m^2.5/mol^0.5/s (LIONSIMBA) to Am^2.5/mol^1.5 (present)
NE.cs_max = 30555
NE.cs0 = 22405
NE.rs = 2e-6
NE.sig = 100
NE.Eac_D = 5000
NE.Eac_k = 5000.
NE.alpha = 0.5
NE.U = x->  0.7222 .+ 0.1387*x + 0.029*x.^0.5 - 0.0172./x + 0.0019./x.^1.5 + 0.2808*exp.(0.9 .- 15*x) - 0.7984*exp.(0.4465*x .- 0.4108)
NE.dUdT = x->  0.001* (0.005269056 .+ 3.299265709*x - 91.79325798*x.^2 + 1004.911008*x.^3 - 5812.278127*x.^4 + 19329.7549*x.^5 - 37147.8947*x.^6 + 38379.18127*x.^7-16515.05308*x.^8) ./
     (1 .- 48.09287227*x + 1017.234804*x.^2 - 10481.80419*x.^3 + 59431.3*x.^4 - 195881.6488*x.^5 + 374577.3152*x.^6 - 385821.1607*x.^7 + 165705.8597*x.^8)

# Electrolyte
EL = Electrolyte()
EL.De = (c_e, T=298)->  10 .^ (-8.43 .- 54 ./ (T .- 229 .- 5e-3 * c_e) - 0.22e-3 * c_e)
EL.kappa = (c_e, T=298)->  1e-4 * c_e .* ( (-10.5 .+ 0.668 * 1e-3 * c_e + 0.494 * 1e-6 * c_e.^2) .+ (0.074 .- 1.78 * 1e-5 * c_e - 8.86 * 1e-10 * c_e.^2) * T .+ (-6.96 * 1e-5 .+ 2.8 * 1e-8 * c_e) * T.^2) .^ 2
EL.dlnf_dlnc = x-> 1
EL.rho = 1290
EL.heat_Q = 134.1
EL.tplus = 0.364
EL.ce0 = 1000

# Separator
SP = Separator()
SP.thickness = 25e-6
SP.lambda = 0.16
SP.rho = 1100
SP.heat_Q = 700
SP.eps = 0.724
SP.eps_fi = 0.
SP.brugg = 4

# Tab
tab = Tab()
tab.width = 40e-3
tab.length = 0.75 * 99.06e-3
tab.area = tab.width * tab.length * 2

# Cell
cell = Cell()    
cell.length = 332.74e-3
cell.width = 99.06e-3
cell.wrapper = 160e-6
cell.I1C = 60
cell.no_layers = 49
cell.capacity = 60
cell.cooling_surface = tab.area
cell.v_h = 4.3
cell.v_l = 2.5
cell.volume = 0.01 * cell.width * cell.length
cell.area =  2.0527
cell.mass = 8.521789028721335e-01
cell.rho = cell.mass / cell.volume  # value for lumped model
cell.alphaT = 0.
cell.heat_Q = 521.7787 # 1464.8 # value for lumped model & check later
cell.h = 150.
cell.T0 = 298
cell.T_amb = cell.T0

# Positive Current Collector
PCC = CurrentCollector()
PCC.thickness = 10e-6
PCC.lambda = 237.
PCC.rho = 2700.
PCC.heat_Q = 897.
PCC.sig =3.55e7

# Negative Current Collector
NCC = CurrentCollector()
NCC.thickness = 10e-6
NCC.lambda = 401.
NCC.rho = 8940.
NCC.heat_Q = 385.
NCC.sig = 5.96e7



# Binder
binder = Binder()
# binder.rho = 

# Scale
scale = Scale()

# assemble to "param_dim"
param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)