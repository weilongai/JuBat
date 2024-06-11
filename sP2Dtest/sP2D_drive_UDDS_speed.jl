using Plots, CSV, DataFrames, Interpolations, BenchmarkTools
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
opt.dt = [1, 20]
opt.dtType = "fixed"
opt.outputType = "manual"
times = zeros(6,3)
volts = zeros(1371,18)

models = ["SPMe", "P2D","sP2D"]
colors = [:black, :blue, :red]
pts=[10, 20, 40, 80, 160, 320]
for j = 1:6
    for k = 1:3
        if k == 1
            opt.Np = 10
            opt.Ns = 10
            opt.Nn = 10
            opt.Nrp = pts[j] 
            opt.Nrn = pts[j] 
        else 
            opt.Np = pts[j]
            opt.Ns = pts[j]
            opt.Nn = pts[j]
            opt.Nrp = 10
            opt.Nrn = 10
        end
        opt.model = models[k]
        case = JuBat.SetCase(param_dim, opt)
        result = JuBat.Solve(case)
        times[j,k] = @belapsed JuBat.Solve($case)
        show(times[j,k])
        plot!(result["time [s]"], result["cell voltage [V]"], label=models[k] * " (JuBat)", linecolor=colors[k])
        volts[:,(k-1) * 6 + j] = result["cell voltage [V]"]
    end
end
path = pwd() * "/src/data/"
result = CSV.read(path * "pybamm_drive_P2D_UDDS.csv", DataFrame, header = 1)
pybamm_1C = Matrix(result)
plot!(pybamm_1C[:,1], pybamm_1C[:,3],label="P2D (PyBaMM)", linestyle =:dot, linecolor=:yellow, lw=2)
plot!(xlabel="time [s]", ylabel="cell voltage [V]", lw=1.5, fontsize=12, size=(400,300), title="1 C discharge", legend=false)
path="c:/code/sP2D/"
savefig(path * "speed_test_UDDS.pdf")
path="c:/code/sP2D/"
#CSV.write(path*"sP2D_UDDS.csv", df, bufsize=10^9)
#savefig("2.png")
# CSV.write("jubat_P2D_$Crate"+"C.csv", result, bufsize=10^9)

# save running time results
times = round.(times, sigdigits = 4)
df =  DataFrame(times,:auto)
headers=["SPMe", "P2D", "sP2D for pts=[10, 20, 40, 80, 160, 320]"]
CSV.write(path * "sP2D_speed_sin_UDDS.csv", df, bufsize=10^9)

volts = round.(volts, sigdigits = 4)
df =  DataFrame(volts,:auto)
headers=["SPMe", "P2D", "sP2D for pts=[10, 20, 40, 80, 160, 320]"]
CSV.write(path * "sP2D_speed_volt_sin_UDDS.csv", df, bufsize=10^9)
