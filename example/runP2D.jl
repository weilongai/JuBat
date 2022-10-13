using Plots, CSV
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
i=-5
opt.Current = x-> i
opt.time = [0 abs(5/i)*3500]
t1 = time()
opt.model = "SPMe"
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)

t2 = time()
opt.model = "SPM"
case1 = JuBat.SetCase(param_dim, opt)
result1 = JuBat.Solve(case1)

t3 = time()
opt.model = "P2D"
case1 = JuBat.SetCase(param_dim, opt)
result2 = JuBat.Solve(case1)
t4 = time()

# Plot
plot(result["time [s]"], result["cell voltage [V]"], label="SPM", lw=1)
plot!(result1["time [s]"], result1["cell voltage [V]"], label="SPMe", lw=1)
plot!(result2["time [s]"], result2["cell voltage [V]"], label="P2D", lw=1)
print("running time: SPM=$(t2-t1) s; SPMe=$(t3-t2) s; P2D=$(t4-t3) s \n")
savefig("2.png")
# CSV.write("result2.csv", result2, bufsize=10^9)