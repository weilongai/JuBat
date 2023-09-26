using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
rates = [0.1, 0.5, 1, 2]
ts = [3700, 3700, 3600, 3600]
opt.Nn = 20
opt.Np = 20
opt.Nrp = 20
opt.Nrn = 20
opt.Ns = 20
for j =1:4
    Crate = rates[j]
    opt.dt = [1, 20]
    opt.dtType = "fixed"
    i= 5*Crate
    opt.Current = x-> i
    opt.time = [ 0 ts[j]/Crate]
    opt.model = "sP2D" # choose model, other options are "SPM" or "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    # recording results
    path="c:/code/sP2D/sP2D_quad" * string(Crate) * 'C'
    times = result["time [s]"]
    voltage = result["cell voltage [V]"]
    phie = result["electrolyte potential [V]"]'
    cel = result["electrolyte lithium concentration [mol/m^3]"]'
    out1 = cat(times, voltage, phie, cel, dims=2)
    out1 = round.(out1, sigdigits = 4)
    df =  DataFrame(out1,:auto)
    headers=["time", "voltage", "phie[:,end]","ce[:,end]"]
    CSV.write(path * ".csv", df, bufsize=10^9)
end