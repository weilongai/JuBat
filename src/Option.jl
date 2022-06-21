using Parameters

@with_kw mutable struct Option
#   option for a lithium-ion battery model
    Np::Int64 = 10
    Ns::Int64 = 10
    Nn::Int64 = 10
    Nrp::Int64 = 10
    Nrn::Int64 = 10
    model::String  = "SPM"
    time::Array{Float64} = [0 3600]
    meshType::String  = "L2"
    gsOrder::Int64 = 2
    dimension::Int64 = 1
    #opt.load = {"constant discharge 1C for 1h"}
    Current::Function = x-> 0
    coupleMethod:: String  = "fully coupled"
    coupleOrder::Int64 = 0
    y0::Array{Float64} = []
    dt::Array{Float64} = [2, 50]
    dtType::String  = "constant" # auto or manual
    dtThreshold::Float64 = 0.01
    solveType::String  = "Crank-Nicolson" # forward, backward or Crank-Nicolson
    outputType::String  = "auto" # auto or manual
    outputTime::Array{Float64,1} = []
    jacobi::String = "constant" # constant or update
end