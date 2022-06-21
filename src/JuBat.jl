module JuBat
using LinearAlgebra, SparseArrays, Plots, Parameters

include("Option.jl") 
include("SetMesh.jl") 
include("SetParams.jl") 
include("SetCase.jl")
include("Assemble.jl") 
include("ElectrodeDiffusion.jl")
include("ElectrolyteDiffusion.jl")
include("SPM.jl") 
include("SPMe.jl") 
include("Solve.jl")
include("Postprocessing.jl")
include("Tools.jl")
include("Thermal.jl")
include("Variables.jl")


export Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell
export Mesh1D, GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables
export SPM, Solve, SPMe
export Arrhenius, IntV
end