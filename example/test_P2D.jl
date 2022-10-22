using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
i=-5
opt.Current = x-> i
opt.time = [0 abs(5/i)*3500]
opt.dtType = "auto"
t3 = time()
opt.model = "P2D"
case1 = JuBat.SetCase(param_dim, opt)
result2 = JuBat.Solve(case1)
t4 = time()

# Plot
plot!(result2["time [s]"], result2["cell voltage [V]"], label="P2D", lw=1)
print("running time:  P2D=$(t4-t3) s \n")
savefig("2.png")
# CSV.write("result2.csv", result2, bufsize=10^9)