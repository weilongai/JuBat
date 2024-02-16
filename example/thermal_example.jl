using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("Northrop")
param_dim.cell.v_l = 2.7
opt = JuBat.Option()
opt.thermalmodel = "lumped"
colors = [:black, :blue, :red, :green]
path = pwd() * "/src/data/"
opt.dtType = "constant"
opt.jacobi = "update"
Crates = [1, 3, 5] #Crates = [1, 3, 5]
I1C = 30 * 2.05272
pV = plot(xlabel="Output capacity [Ah]", ylabel="Cell voltage [V]", legend_column = 2)
pT = plot(xlabel="Output capacity [Ah]", ylabel="Temperature [K]", legend_column = 2)
for i in eachindex(Crates)
    opt.Current = x->  I1C * Crates[i]
    if i < 3
        opt.dt = [0.1, 0.1]
    else
        opt.dt = [0.01, 0.01]
    end
    opt.time = [0 3540/Crates[i]]
    opt.model = "SPMe" # choose model, other options are "SPM" or "SPMe"
    case1 = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case1)
    plot!(pV, result["time [s]"]/3600*Crates[i]*I1C, result["cell voltage [V]"], label=string(Crates[i]) * "C (JuBat)", linecolor=colors[i], legend_column = 2)
    plot!(pT, result["time [s]"]/3600*Crates[i]*I1C, result["temperature [K]"], label=string(Crates[i]) * "C (JuBat)", linecolor=colors[i], legend_column = 2)
    result1 = CSV.read(path * "lionsimba" * string(Crates[i]) * "c.csv", DataFrame, header = 1)
    lionsimba = Matrix(result1)
    plot!(pV, lionsimba[:,1]/3600*Crates[i]*I1C, lionsimba[:,3],label=string(Crates[i]) * "C (LIONSIMBA)", linestyle =:dot, linecolor=colors[i], lw=2, legend_column = 2)
    plot!(pT, lionsimba[:,1]/3600*Crates[i]*I1C, lionsimba[:,4],label=string(Crates[i]) * "C (LIONSIMBA)", linestyle =:dot, linecolor=colors[i], lw=2, legend_column = 2)
end

savefig(pV, "thermal example-V.pdf")
savefig(pT, "thermal example-T.pdf")
