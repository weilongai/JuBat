using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
opt.thermalmodel = "lumped"
Crate = 2          
I1C = 5
i= I1C * Crate
opt.Current = x-> i
opt.dtType = "constant"
opt.dt = [0.2, 0.2]
opt.time = [0 3500]
opt.model = "P2D" # choose model, other options are "SPM" or "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)
plot(result["time [s]"], result["cell voltage [V]"], label="P2D", xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
savefig("minimal_example_test1.pdf")
plot(result["time [s]"], result["temperature [K]"], label="P2D", xlabel="time [s]", ylabel="Temperature [K]", lw=1)
savefig("minimal_example_test2.pdf")