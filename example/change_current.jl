using Plots, CSV, DataFrames, Interpolations
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
path = pwd()
file= path * "/src/data/drive_cycles/US06.csv"
data = CSV.read(file, DataFrame, header = 2)
current=Matrix(data)
current_intpolation = LinearInterpolation(current[:,1],current[:,2])
opt.Current = x-> current_intpolation(x)
opt.time = collect(0:1:600)
opt.dt = [0.1, 12]
opt.model = "P2D"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)

# Plot
times = result["time [s]"]
voltage = result["cell voltage [V]"]
current = result["cell current [A]"]
plot(times, voltage, label="P2D",xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
