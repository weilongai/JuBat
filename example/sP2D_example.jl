using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
rates = [ 0.5, 1, 2]
ts = [3700, 3700, 3600, 3600]
for j =1:3
    Crate = rates[j]
    opt.dt = [1, 20]/Crate
    opt.dtType = "fixed"
    i= 5*Crate
    opt.Current = x-> i
    opt.time = [ 0 ts[j]/Crate]
    opt.model = "sP2D" # choose model, other options are "SPM" or "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    opt.model = "P2D"
    case2 = JuBat.SetCase(param_dim, opt)
    result2 = JuBat.Solve(case2)

    # ploting results
    plot!(result["time [s]"], result["cell voltage [V]"], label="sP2D", xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
    plot!(result2["time [s]"], result2["cell voltage [V]"], label="P2D", linestyle =:dot, linecolor =:black, lw=2.5)
end
savefig("sP2D_validation.pdf")
JuBat.Citation(["ai2024b"])