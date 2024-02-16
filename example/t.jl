using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_l = 2.7
opt = JuBat.Option()
opt.thermalmodel = "lumped"
colors = [:black, :blue, :red, :green]
path = pwd() * "/src/data/"
opt.dtType = "constant"
opt.jacobi = "update"
Crates = [1, 2] #Crates = [1, 3, 5]
I1C = 5
pV = plot(xlabel="Output capacity [Ah]", ylabel="Cell voltage [V]", legend_column = 2)
pT = plot(xlabel="Output capacity [Ah]", ylabel="Temperature [K]", legend_column = 2)
for i in eachindex(Crates)
    result1 = CSV.read(path * "pybamm_P2D_" * string(Crates[i]) * "C_T.csv", DataFrame, header = 1)
    result1 = Matrix(result1)
    plot!(pV, result1[:,1]/3600*Crates[i]*30, result1[:,3],label=string(Crates[i]) * "C (PyBaMM)", linestyle =:dot, linecolor=colors[i], lw=2, legend_column = 2)
    plot!(pT, result1[:,1]/3600*Crates[i]*30, result1[:,4],label=string(Crates[i]) * "C (PyBaMM)", linestyle =:dot, linecolor=colors[i], lw=2, legend_column = 2)
end

savefig(pV, "thermal example-V.pdf")
savefig(pT, "thermal example-T.pdf")
