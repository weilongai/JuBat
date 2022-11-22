using Plots
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
i= 5
opt.Current = x-> i
opt.time = [0 abs(5/i)*3600]
t1 = time()
case = JuBat.SetCase(param_dim, opt)
result = JuBat.SolveODE(case)
t2 = time()
# opt.model = "SPMe"
# case1 = JuBat.SetCase(param_dim, opt)
# result1 = JuBat.Solve(case1)
t3 = time()
# opt.model = "P2D"
# case1 = JuBat.SetCase(param_dim, opt)
# result2 = JuBat.Solve(case1)

# Plot
plot(result["time [s]"], result["cell voltage [V]"])
# plot!(result1["time [s]"], result1["cell voltage [V]"])
# plot!(result2["time [s]"], result2["cell voltage [V]"])
print("running time: SPM=$(t2-t1) s; SPMe=$(t3-t2) s; \n")
savefig("1.png")