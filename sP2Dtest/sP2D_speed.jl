using Plots, CSV, DataFrames, BenchmarkTools
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
i = 5
opt.Current = x-> i
opt.dt = [1, 1]
opt.time = [0 abs(5/i)*3540]
pts=[10, 20, 40, 80, 160, 320]
models = ["SPMe", "P2D","sP2D"]
colors = [:black, :blue, :red]
times = zeros(6,3)
volts = zeros(3541,18)
path = pwd() * "/src/data/"
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
        # @belapsed JuBat.Solve($case)
        # bt = @benchmark JuBat.Solve($case)
        # times[j,k] = mean(bt.times)
        # show(bt)
        times[j,k] = @belapsed JuBat.Solve($case)
        show(times[j,k])
        plot!(result["time [s]"], result["cell voltage [V]"], label=models[k] * " (JuBat)", linecolor=colors[k])
        volts[:,(k-1) * 6 + j] = result["cell voltage [V]"]
    end
end
result = CSV.read(path * "pybamm_P2D_1C.csv", DataFrame, header = 0)
pybamm_1C = Matrix(result)
plot!(pybamm_1C[:,1], pybamm_1C[:,2],label="P2D (PyBaMM)", linestyle =:dot, linecolor=:yellow, lw=2)
plot!(xlabel="time [s]", ylabel="cell voltage [V]", lw=1.5, fontsize=12, size=(400,300), title="1 C discharge", legend=:bottomleft)
path="c:/code/sP2D/"
savefig(path * "speed_test.pdf")

# save running time results
times = round.(times, sigdigits = 4)
df =  DataFrame(times,:auto)
headers=["SPMe", "P2D", "sP2D for pts=[10, 20, 40, 80, 160, 320]"]
CSV.write(path * "sP2D_speed_quad.csv", df, bufsize=10^9)

volts = round.(volts, sigdigits = 4)
df =  DataFrame(volts,:auto)
headers=["SPMe", "P2D", "sP2D for pts=[10, 20, 40, 80, 160, 320]"]
CSV.write(path * "sP2D_speed_volt_quad.csv", df, bufsize=10^9)

 