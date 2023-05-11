using Plots, CSV, DataFrames, Interpolations
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
path =pwd()
file= path * "/src/data/drive_cycles/UDDS.csv"
data = CSV.read(file, DataFrame, header = 3)
current=Matrix(data)
current_intpolation = LinearInterpolation(current[:,1],current[:,2])
opt.Current = x-> current_intpolation(x)
opt.time = collect(0:1:1369)
opt.dt = [0.1, 12]
opt.dtType = "auto"
opt.outputType = "manual"

opt.model = "sP2D"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)

# Plot
plot!(result["time [s]"], result["cell voltage [V]"], label="P2D", lw=1)
times = result["time [s]"]
voltage = result["cell voltage [V]"]
current = result["cell current [A]"]
out = cat(times, current, voltage,dims=2)
df =  DataFrame(out,:auto)
path="c:/code/sP2D/"
CSV.write(path*"sP2D_UDDS.csv", df, bufsize=10^9)
#savefig("2.png")
# CSV.write("jubat_P2D_$Crate"+"C.csv", result, bufsize=10^9)