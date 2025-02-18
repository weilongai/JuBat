"""
Parameters for the Enertech cell (Ai2020)
"""
# Positive Electrode
PE = Electrode()
PE.theta_100 = 0.435 #
PE.theta_0 = 0.9651 #
PE.thickness = 6.8e-5 # 
PE.lambda = 1.58 #
PE.Ds = 5.387e-15 #
PE.rho = 2470 #
PE.heat_Q = 1080.2 #
PE.eps =  0.32 #
PE.eps_fi = 0
PE.brugg = 1.83 #
PE.k = 1e-11 * 96485.33289 #
PE.cs_max = 49943 #
PE.cs0 = 21725 #
PE.rs = 3e-6 #
PE.sig = 10. #
PE.Eac_D = 5000 #
PE.Eac_k = 5000 #
PE.alpha = 0.5 #

PE.U = x->(-107897.40 *x.^9  
.+677406.28 *x.^8
.-1873803.91 *x.^7
.+2996535.44 *x.^6
.-3052331.36 *x.^5
.+2053377.31 *x.^4
.-912135.88 *x.^3
.+257964.35 *x.^2
.-42146.98 *x
.+3035.67)
PE.dUdT = x-> (
   .-3.20392657 * x.^7
      .+ 14.5719049 * x.^6
      .-27.9047599 * x.^5
      .+ 29.1744564 * x.^4
      .-17.992018 * x.^3
      .+ 6.54799331 * x.^2
      .-1.30382445* x
      .+ 0.109667298
  )

# Negative Electrode
NE = Electrode()
NE.theta_100 = 0.84 #
NE.theta_0 = 0.0065 #
NE.thickness = 7.65e-5 # 
NE.lambda = 1.04 #
NE.Ds = 3.9e-14 #
NE.rho = 2470 #
NE.heat_Q = 1080.2 #
NE.eps = 0.33
NE.eps_fi = 0
NE.brugg = 2.914 #
NE.k = 1e-11 * 96485.33289 #
NE.cs_max = 29700 #
NE.cs0 = 24108 #
NE.rs = 5e-6
NE.sig = 100. #
NE.Eac_D = 5000. #
NE.Eac_k = 5000. #
NE.alpha = 0.5 #
NE.U = x-> (-2058.29865*x.^9
.+ 10040.08960*x.^8
.-20824.86740*x.^7
.+23911.86578*x.^6
.-16576.3692*x.^5
.+7098.09151*x.^4
.-1845.43634*x.^3
.+275.31114*x.^2
.-21.20097*x
.+ 0.84498)
NE.dUdT = x->(
   0.001
   * (
       0.005269056
       .+ 3.299265709 * x
       .- 91.79325798 * x.^2
       .+ 1004.911008 * x.^3
       .- 5812.278127 * x.^4
       .+ 19329.7549 * x.^5
       .- 37147.8947 * x.^6
       .+ 38379.18127 * x.^7
       .- 16515.05308 * x.^8
   )
   ./ (
       1
       .- 48.09287227 * x
       .+ 1017.234804 * x.^2
       .- 10481.80419 * x.^3
       .+ 59431.3 * x.^4
       .- 195881.6488 * x.^5
       .+ 374577.3152 * x.^6
       .- 385821.1607 * x.^7
       .+ 165705.8597 * x.^8
   )
)


# Electrolyte
EL = Electrolyte()
EL.De = (x, y=0)-> 10 .^ (-4.43 .- 54 ./ (298.15 - 229 .- 5e-3 * x) .- 0.22e-3 * x)
EL.kappa = (x, y=0)-> 1e-4* (x).* (
   (-10.5 .+ 0.668 * 1e-3 .* (x) .+ 0.494 * 1e-6 .*(x).^2)
   .+ (0.074 .- 1.78 * 1e-5 .* (x) .- 8.86 * 1e-10 .* (x).^2) * 298.15
   .+ (-6.96 * 1e-5 .+ 2.8 * 1e-8 .* (x)) * 298.15.^2
).^ 2
EL.dlnf_dlnc = x-> 1
EL.rho = 1290 #？
EL.heat_Q = 0 #
EL.tplus = 0.38 #
EL.ce0 = 1000 #

# Separator
SP = Separator()
SP.thickness = 2.5e-5 #
SP.lambda = 0.334 #
SP.rho = 2470 #
SP.heat_Q = 1080.2 #
SP.eps = 0.5 #
SP.eps_fi = 0.
SP.brugg = 1.5  #

# Cell
cell = Cell()    
cell.length = 0.051 #
cell.width = 0.047 #
# cell.wrapper = 0.
cell.I1C = 2.28 #
cell.no_layers = 34 #
cell.capacity = 2.28 #
cell.cooling_surface = 0.0060484 #
cell.v_h = 4.2 #
cell.v_l = 3.0 #
cell.volume = 1.5341e-5 #
cell.rho = 0.
cell.alphaT = 0.
cell.heat_Q = 0.
cell.h = 0.
cell.T0 = 298.15 #

# Positive Current Collector
PCC = CurrentCollector()
PCC.thickness = 1.5e-5 #
PCC.lambda = 237. #
PCC.rho = 2700. #
PCC.heat_Q = 897. #
PCC.sig =3.6914e7 #

# Negative Current Collector
NCC = CurrentCollector()
NCC.thickness = 1e-5 #
NCC.lambda = 401. #
NCC.rho = 8960. #
NCC.heat_Q = 385. #
NCC.sig = 5.8411e7 #

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