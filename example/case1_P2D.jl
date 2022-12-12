using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
Crate = 1
i= 5 * Crate
opt.Current = x-> i
opt.time = collect(0:10:3600)
opt.dt = [0.1, 12]
opt.dtType = "auto"
opt.outputType = "manual"
opt.gsorder = 2
Nx = 20
opt.Np = Nx
opt.Ns = Nx
opt.Nn = Nx
opt.Nrp = Nx
opt.Nrn = Nx
t3 = time()
opt.model = "P2D"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)
t4 = time()

# Plot
plot!(result["time [s]"], result["cell voltage [V]"], label="P2D", lw=1)
print("running time:  P2D=$(t4-t3) s \n")
times = result["time [s]"]
csn = result["negative particle surface lithium concentration [mol/m^3]"]'
csp = result["positive particle surface lithium concentration [mol/m^3]"]'
cel = result["electrolyte lithium concentration [mol/m^3]"]'
phiel = result["electrolyte potentials [V]"]'
v=[1, 101, 201, 355]
out = cat(times[v], csn[v,:], csp[v,:], cel[v,:], phiel[v,:],dims=2)
df =  DataFrame(out,:auto)
path="c:/code/"
CSV.write(path*"jubat_P2D_$(Crate)C_details.csv", df, bufsize=10^9)
#savefig("2.png")
# CSV.write("jubat_P2D_$Crate"+"C.csv", result, bufsize=10^9)
