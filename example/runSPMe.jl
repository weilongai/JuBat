using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
Crate = 2
i= 5 * Crate
opt.Current = x-> i
opt.time = [0 3600/Crate]
opt.dt = [2/Crate, 4/Crate]
opt.dtType = "auto"
t1 = time()
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)
t2 = time()
opt.model = "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result1 = JuBat.Solve(case1)
t3 = time()
# opt.model = "P2D"
# case1 = JuBat.SetCase(param_dim, opt)
# result2 = JuBat.Solve(case1)

# Plot
plot(result["time [s]"], result["cell voltage [V]"], label="SPM")
plot!(result1["time [s]"], result1["cell voltage [V]"], label="SPMe")
# plot!(result2["time [s]"], result2["cell voltage [V]"])
print("running time: SPM=$(t2-t1) s; SPMe=$(t3-t2) s; \n")
savefig("1.png")
path="c:/code/"
out = cat(result["time [s]"], result["cell voltage [V]"],dims=2)
df =  DataFrame(out,:auto)
CSV.write(path * "jubat_SPM_$(Crate)C.csv", df, bufsize=10^9)
out1 = cat(result1["time [s]"], result1["cell voltage [V]"],dims=2)
df =  DataFrame(out1,:auto)
CSV.write(path * "jubat_SPMe_$(Crate)C.csv", df, bufsize=10^9)