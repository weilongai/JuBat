using Plots, CSV
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
i = 5
opt.Current = x-> i
opt.time = collect(0:100:abs(5/i)*3500)
t1 = time()
opt.model = "SPM"
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)

t2 = time()
opt.dtType = "auto"
case1 = JuBat.SetCase(param_dim, opt)
result1 = JuBat.Solve(case1)
t3 = time()

# Plot
plot(result["time [s]"], result["cell voltage [V]"], label="constant dt", xlabel="time [s]", ylabel="cell voltage [V]",lw=1)
plot!(result1["time [s]"], result1["cell voltage [V]"], label="auto dt", lw=1) # shape=:circle
print("running time: constant dt=$(t2-t1) s; auto dt=$(t3-t2) s; \n")

# CSV.write("result2.csv", result2, bufsize=10^9)