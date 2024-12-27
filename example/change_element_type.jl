using Plots, CSV
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
i = 5
opt.Current = x-> i
opt.time = collect(0:100:abs(5/i)*3500)
opt.model = "P2D"
opt.meshType  = "L2"
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)

opt.meshType  = "L3"
case1 = JuBat.SetCase(param_dim, opt)
result1 = JuBat.Solve(case1)

# Plot
plot(result["time [s]"], result["cell voltage [V]"], label="Linear element", xlabel="time [s]", ylabel="cell voltage [V]",lw=1.5)
plot!(result1["time [s]"], result1["cell voltage [V]"], label="Quadratic element", linestyle=:dot, lw=2) # shape=:circle
savefig("change_element.pdf")
JuBat.Citation(["ai2023"])
# CSV.write("result2.csv", result2, bufsize=10^9)