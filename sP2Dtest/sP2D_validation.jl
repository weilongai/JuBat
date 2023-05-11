using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
opt.dt = [1, 20]
opt.dtType = "fixed"
Crate = 2
i= 5*Crate
opt.Current = x-> i
opt.time = [ 0 3500/Crate]
opt.model = "sP2D" # choose model, other options are "SPM" or "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)

opt.model = "P2D"
case2 = JuBat.SetCase(param_dim, opt)
result2 = JuBat.Solve(case2)

# ploting results
plot(result["time [s]"], result["cell voltage [V]"], label="sP2D", xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
plot!(result2["time [s]"], result2["cell voltage [V]"], label="P2D", linestyle =:dot, linecolor =:black, lw=2.5)
savefig("sP2D.png")

mesh_el = case1.mesh["electrolyte"]
v=Int(round(opt.time[2]/5))*[1, 2, 3, 4, 5]
plot()
for i in eachindex(v)
    plot!(mesh_el.node, result["electrolyte potential [V]"][:,v[i]], label="sP2D", xlabel="x/L", ylabel="potential [V]", lw=1)
    plot!(mesh_el.node, result2["electrolyte potential [V]"][:,v[i]], label="P2D", linestyle =:dot, linecolor =:black, lw=2.5)
end 
savefig("sP2D_phie.png")

plot()
for i in eachindex(v)
    plot!(mesh_el.node, result["electrolyte potential [V]"][:,v[i]] .- result["electrolyte potential [V]"][16,v[i]], label="sP2D", xlabel="x/L", ylabel="potential [V]", lw=1)
    plot!(mesh_el.node, result2["electrolyte potential [V]"][:,v[i]] .- result2["electrolyte potential [V]"][16,v[i]], label="P2D", linestyle =:dot, linecolor =:black, lw=2.5)
end 
savefig("sP2D_phie_rel.png")

plot()
for i in eachindex(v)
    plot!(mesh_el.node, result["electrolyte lithium concentration [mol/m^3]"][:,v[i]], label="sP2D", xlabel="x/L", ylabel="lithium concentration [mol/m^3]", lw=1)
    plot!(mesh_el.node, result2["electrolyte lithium concentration [mol/m^3]"][:,v[i]], label="P2D", linestyle =:dot, linecolor =:black, lw=2.5)
end 
savefig("sP2D_phie_ce.png")