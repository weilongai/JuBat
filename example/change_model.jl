using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
i = 5
opt.Current = x-> i
opt.time = [0 abs(5/i)*3540]
models = ["SPM", "SPMe", "P2D"]
colors = [:black, :blue, :red]
path = pwd() * "/src/data/"
for i = 1:3
    opt.model = models[i]
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
    plot!(result["time [s]"], result["cell voltage [V]"], label=models[i] * " (JuBat)", linecolor=colors[i])
    result = CSV.read(path * "pybamm_" * models[i] * "_1C.csv", DataFrame, header = 0)
    pybamm_1C = Matrix(result)
    plot!(pybamm_1C[:,1], pybamm_1C[:,2],label=models[i] * " (PyBaMM)", linestyle =:dot, linecolor=colors[i], lw=2)
end
plot!(xlabel="time [s]", ylabel="cell voltage [V]", lw=1.5, fontsize=12, size=(400,300), title="1 C discharge", legend=:bottomleft)
savefig("change_model.pdf")