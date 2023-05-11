using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
Crate = 0.1
opt.dt = [1, 20]
opt.dtType = "fixed"
i= 5*Crate
opt.Current = x-> i
opt.time = [ 0 3700/Crate]
opt.model = "sP2D" # choose model, other options are "SPM" or "SPMe"
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)

# recording results
path="c:/code/sP2D/sP2D_quad" * string(Crate) * 'C'
times = result["time [s]"]
voltage = result["cell voltage [V]"]
phie = result["electrolyte potential [V]"]'
cel = result["electrolyte lithium concentration [mol/m^3]"]'
out1 = cat(times, voltage, phie, cel, dims=2)
out1 = round.(out1, sigdigits = 4)
df =  DataFrame(out1,:auto)
headers=["time", "voltage", "phie[:,end]","ce[:,end]"]
CSV.write(path * ".csv", df, bufsize=10^9)