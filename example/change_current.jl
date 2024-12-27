using Plots, CSV, DataFrames, Interpolations
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_h = 4.3
opt = JuBat.Option()
path = pwd() * "/src/data/drive_cycles"
data = CSV.read(path * "/UDDS.csv", DataFrame, header = 3)
current=Matrix(data)
current_intpolation = LinearInterpolation(current[:,1],current[:,2])
opt.Current = x-> current_intpolation(x)
opt.time = collect(0:1:1369)
opt.dt = [0.5, 12]
opt.model = "P2D"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)
times = result["time [s]"]
voltage = result["cell voltage [V]"]
current = result["cell current [A]"]
plot(times, voltage, label="JuBat",xlabel="time [s]", ylabel="cell voltage [V]",linecolor=:black, lw=1.5, fontsize=12, size=(400,300))
path = pwd() * "/src/data/"
result = CSV.read(path * "pybamm_drive_P2D_UDDS.csv", DataFrame, header = 1)
pybamm_drive = Matrix(result)
plot!(pybamm_drive[:,1], pybamm_drive[:,3],label="PyBaMM", linestyle=:dot, linecolor=:red, lw=2, title="Drive cycle UDDS", legend=(0.8,0.95))
savefig("drive_UDDS_V.pdf")
JuBat.Citation()
# plot(times, current,xlabel="time [s]", ylabel="current [A]", label="Current", fontsize=12, size=(400,300), title="Drive cycle UDDS")
# savefig("drive_UDDS_I.pdf")