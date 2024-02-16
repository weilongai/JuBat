using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("Northrop")
opt = JuBat.Option()
opt.thermalmodel = "none"
Crate = 0.1
I1C = 30 * 2.05272
i= I1C * Crate
opt.Current = x-> i
opt.dtType = "constant"
opt.jacobi = "update"
opt.dt = [0.2, 0.2]
opt.time = [0 3500]
opt.model = "P2D" # choose model, other options are "SPM" or "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)
plot(result["time [s]"], result["cell voltage [V]"], label="P2D", xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
savefig("thermal example-V.pdf")
plot(result["time [s]"], result["temperature [K]"], label="P2D", xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
savefig("thermal example-T.pdf")
