using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
Crate = 2
opt.dt = [1, 10]
opt.dtType = "fixed"
i= 5*Crate
opt.Nn = 20
opt.Np = 20
opt.Nrp = 20
opt.Nrn = 20
opt.Ns = 20
opt.Current = x-> i
opt.time = [ 0 3500/Crate]
opt.model = "P2D" # choose model, other options are "SPM" or "SPMe"
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)

param_dim1 = deepcopy(param_dim)
param_dim2 = deepcopy(param_dim)
param_dim1.EL.kappa = x-> param_dim.EL.kappa(x) / 2
param_dim2.EL.kappa = x-> param_dim.EL.kappa(x) * 2
case1 = JuBat.SetCase(param_dim1, opt)
result1 = JuBat.Solve(case1)
case2 = JuBat.SetCase(param_dim2, opt)
result2 = JuBat.Solve(case2)

# recording results
path="c:/code/sP2D/prestudy/P2D_" * string(Crate) * "C_kappa"
times = result["time [s]"]
voltage = result["cell voltage [V]"]
phie = result["electrolyte potential [V]"]'
cel = result["electrolyte lithium concentration [mol/m^3]"]'
times1 = result1["time [s]"]
voltage1 = result1["cell voltage [V]"]
phie1 = result1["electrolyte potential [V]"]'
cel1 = result1["electrolyte lithium concentration [mol/m^3]"]'
times2 = result2["time [s]"]
voltage2 = result2["cell voltage [V]"]
phie2 = result2["electrolyte potential [V]"]'
cel2 = result2["electrolyte lithium concentration [mol/m^3]"]'
out1 = cat(times, voltage, phie, cel, times1, voltage1, phie1, cel1, times2, voltage2, phie2, cel2, dims=2)
out1 = round.(out1, sigdigits = 5)
df =  DataFrame(out1,:auto)
headers=["kappa:time", "voltage", "phie[:,end]","ce[:,end]", "kappa/10:time", "voltage", "phie[:,end]","ce[:,end]", "kappa*10:time", "voltage", "phie[:,end]","ce[:,end]"]
CSV.write(path * ".csv", df, bufsize=10^9)