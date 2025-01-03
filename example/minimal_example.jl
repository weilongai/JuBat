using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("Enertech")
opt = JuBat.Option()
Crate = 1
i= 5*Crate
opt.Current = x-> i
opt.time = [0 3600]
opt.model = "P2D" # choose model, other options are "SPM" or "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)
plot(result["time [s]"], result["cell voltage [V]"], label="P2D", xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
JuBat.Citation()
savefig("minimal_example.pdf")
