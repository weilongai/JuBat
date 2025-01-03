using Plots, CSV, DataFrames, Revise
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("Enertech")
opt = JuBat.Option()
Crate = 1
i= 2.28*Crate
opt.Current = x-> i
opt.time = [0 9000]
opt.dtType = "auto"
opt.dtThreshold = 0.0001
opt.model = "P2D" # choose model, other options are "SPM" or "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)

c_vector = vec(result["electrolyte lithium concentration [mol/m^3]"])
c_df = DataFrame(c = c_vector)
file_path = ("C:/Users/user/Desktop/resultenertech/ceori.csv")
CSV.write(file_path, c_df)