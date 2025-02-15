module JuBat
using LinearAlgebra, SparseArrays, Plots, Parameters, CSV, Infiltrator

include("Option.jl") 
include("SetMesh.jl") 
include("SetParams.jl") 
include("SetCase.jl")
include("Assemble.jl") 
include("ElectrodeDiffusion.jl")
include("ElectrolyteDiffusion.jl")
include("ElectrodePotential.jl")
include("ElectrolytePotential.jl")
include("SPM.jl") 
include("SPMe.jl") 
include("P2D.jl") 
include("Solve.jl")
include("PostProcessing.jl")
include("Tools.jl")
include("Thermal.jl")
include("Variables.jl")
include("Initialisation.jl")
include("sP2D.jl")
include("Citation.jl")
include("Mechanical.jl")


export Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell
export Mesh1D, GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables
export SPM, Solve, SPMe
export Arrhenius, IntV
end