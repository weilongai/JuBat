using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
Crate = 1
i= 5*Crate
opt.Current = x-> i
opt.time = [0 abs(5/i)*3500]
opt.dt = [1 1]
opt.dtType = "constant"
opt.model = "SPM"
case = JuBat.SetCase(param_dim, opt)
t1 = time()
result = JuBat.Solve(case)

t2 = time()

# Plot
plot(result["time [s]"], result["cell voltage [V]"], label="jubat_SPM", lw=1)
path =pwd()
file = path * "/example/pybamm_SPM_1C.txt"
result = CSV.read(file, DataFrame, header = 0)
pybamm_SPM_1C = Matrix(result)
plot!(pybamm_SPM_1C[1,:], pybamm_SPM_1C[2,:], label="pybamm_SPM", lw=1.5)
savefig("SPM.png")
# CSV.write("result2.csv", result2, bufsize=10^9)
