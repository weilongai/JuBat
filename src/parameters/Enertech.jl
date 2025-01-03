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
PE.eps_fi = 0.06
PE.brugg = 1.83 #
PE.k = 1e-11 * 96485.33289 #
PE.cs_max = 49943 #
PE.cs0 = 21725 #
PE.rs = 3e-6 #
PE.sig = 10. #
PE.Eac_D = 5000 #
PE.Eac_k = 5000 #
PE.alpha = 0.5 #
PE.E = 375000000000.0
PE.nu = 0.2
PE.Omega = -7.28e-7
PE.U = x->(-107897.40 *x.^9  .+677406.28 *x.^8 .-1873803.91 *x.^7 .+2996535.44 *x.^6 .-3052331.36 *x.^5 .+2053377.31 *x.^4 .-912135.88 *x.^3 .+257964.35 *x.^2 .-42146.98 *x .+3035.67)

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
NE.eps_fi = 0.06
NE.brugg = 2.914 #
NE.k = 1e-11 * 96485.33289 #
NE.cs_max = 28700 #
NE.cs0 = 24108 #
NE.rs = 5e-6
NE.sig = 100. #
NE.Eac_D = 5000. #
NE.Eac_k = 5000. #
NE.alpha = 0.5 #
NE.E = 15000000000.0
NE.nu = 0.3
NE.Omega = 3.1e-6
NE.U = x -> (
   -21333233060.3788 .* x.^25 .+ 
   214414158522.696 .* x.^24 .- 
   967791558850.832 .* x.^23 .+ 
   2551110424909.95 .* x.^22 .- 
   4216223176986.89 .* x.^21 .+ 
   4154322237396.47 .* x.^20 .- 
   1441975443918.20 .* x.^19 .- 
   2145633957566.25 .* x.^18 .+ 
   3558015865082.06 .* x.^17 .- 
   1806361400171.50 .* x.^16 .- 
   1178075855900.44 .* x.^15 .+ 
   2977112532980.41 .* x.^14 .- 
   2925950140441.38 .* x.^13 .+ 
   1901113351779.06 .* x.^12 .- 
   907042269942.721 .* x.^11 .+ 
   329682508054.031 .* x.^10 .- 
   92457648628.5046 .* x.^9 .+ 
   20007802947.4577 .* x.^8 .- 
   3310569828.43730 .* x.^7 .+ 
   411459245.533657 .* x.^6 .- 
   37373632.1267117 .* x.^5 .+ 
   2387565.33605553 .* x.^4 .- 
   102350.934580862 .* x.^3 .+ 
   2849.28977163607 .* x.^2 .- 
   54.6653601804887 .* x .+ 
   0.974968767152609
)
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
EL.De = (x, y=0)-> 10 .^ (-4.43 .- 54 ./ (y - 229 .- 5e-3 * x) .- 0.22e-3 * x)
EL.kappa = (x, y=0)-> 1e-4* (x).* (
   (-10.5 .+ 0.668 * 1e-3 .* (x) .+ 0.494 * 1e-6 .*(x).^2)
   .+ (0.074 .- 1.78 * 1e-5 .* (x) .- 8.86 * 1e-10 .* (x).^2) * 298.15
   .+ (-6.96 * 1e-5 .+ 2.8 * 1e-8 .* (x)) * 298.15^2
).^ 2
EL.dlnf_dlnc = x-> ( 0.601 .- 0.24 * (x./1000).^ 0.5 .+0.982)/(1 - 0.38)
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
cell.area = cell.width * cell.length * cell.no_layers
cell.mass = cell.rho * cell.volume

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
