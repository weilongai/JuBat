using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
path = "C:\\Users\\user\\Desktop\\JuBat\\src\\data\\"
Crate = 1
i= 5*Crate
opt.Current = x-> i
opt.time = [0 3600]
opt.model = "SPM" # choose model, other options are "SPM" or "SPMe"

opt.mechanicalmodel = "full"
case1 = JuBat.SetCase(param_dim, opt)
result1 = JuBat.Solve(case1)
result2 = CSV.read(path * "pybamm_SPM_1C_Stress_t.csv", DataFrame, header = 1)
pybamm_drive = Matrix(result2)

plot(result1["time [s]"], result1["negative particle surface tangential stress[Pa]"]', label="Jubat", xlabel="time [s]", ylabel="negative particle surface tangential stress[Pa]", lw=1)
plot!(pybamm_drive[:,1], pybamm_drive[:,2],label="PyBaMM", linestyle=:dot, linecolor=:red, lw=2, title="mechanical validation", legend=(0.8,0.95))

JuBat.Citation()
savefig("mechanical_example_SPM.pdf")