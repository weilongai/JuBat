module JuBat
using LinearAlgebra, SparseArrays, Plots, Parameters

include("Option.jl") 
include("SetMesh.jl") 
include("SetParams.jl") 
include("SetCase.jl")
include("Assemble.jl") 
include("ElectrodeDiffusion.jl") 
include("SPM.jl") 
include("Solve.jl")
include("Postprocessing.jl")

export Assemble, ElectrodeDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell
export Mesh1D, GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam
export SPM, Solve
end