using CSV, DataFrames, Plots
path =pwd()
file = path * "/example/pybamm_DFN_1C.txt"
result = CSV.read(file, DataFrame, header = 0)
pybamm_DFN_1C = Matrix(result)
plot(pybamm_DFN_1C[1,:], pybamm_DFN_1C[2,:], label="pybamm", lw=1.5)