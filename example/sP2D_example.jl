using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
Crate = 2
i= 5*Crate
opt.Current = x-> i
opt.time = [ 0 3600/Crate]
opt.model = "sP2D" # choose model, other options are "SPM" or "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)

opt.model = "P2D"
case2 = JuBat.SetCase(param_dim, opt)
result2 = JuBat.Solve(case2)
plot(result["time [s]"], result["cell voltage [V]"], label="sP2D", xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
plot!(result2["time [s]"], result2["cell voltage [V]"], label="P2D", linestyle =:dot, linecolor =:black, lw=2.5)
savefig("sP2D.png")