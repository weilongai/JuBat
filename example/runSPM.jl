include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
opt.Current = x-> -5
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)
# Plot