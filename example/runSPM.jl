using Plots
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
i=-25
opt.Current = x-> i
opt.time = [0 abs(5/i)*3600*0.7]
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)

opt.model = "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result1 = JuBat.Solve(case1)


# Plot
plot(result["time [s]"], result["cell voltage [V]"])
plot!(result1["time [s]"], result1["cell voltage [V]"])
savefig("1.png")