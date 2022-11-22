using Plots, CSV, DataFrames, BenchmarkTools
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
Crate = 2
i= 5 * Crate
opt.Current = x-> i
opt.time = [0 900]
opt.dt = [0.1, 0.1]
opt.dtType = "constant"
opt.gsorder = 2
Nx = 10
opt.Np = Nx
opt.Ns = Nx
opt.Nn = Nx
opt.Nrp = Nx
opt.Nrn = Nx
opt.meshType  = "L2"
opt.model = "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
# @btime JuBat.Solve(case1)
t2 = time()
result1 = JuBat.Solve(case1)
t3 = time()


# Plot
plot!(result1["time [s]"], result1["cell voltage [V]"], label="SPMe")
# plot!(result2["time [s]"], result2["cell voltage [V]"])
print("running time: SPMe=$(t3-t2) s; \n")
savefig("1.png")
path="c:/code/jubat_"
times = result1["time [s]"][end]
voltage = result1["cell voltage [V]"][end]
csn = result1["negative particle lithium concentration [mol/m^3]"][:, end]
csp = result1["positive particle lithium concentration [mol/m^3]"][:, end]
cel = result1["electrolyte lithium concentration [mol/m^3]"][:, end]
out1 = vcat(times, voltage, csn, csp, cel)
df =  DataFrame([out1],:auto)
headers=["time", "voltage", "csn[:,end]","csp[:,end]", "cel[:,end]"]
CSV.write(path * "Nx$(Nx)_SPMe_2C.csv", df, bufsize=10^9, header = headers)