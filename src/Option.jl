using Parameters

@with_kw mutable struct Option
#   option for a lithium-ion battery model
    Np::Int32 = 10
    Ns::Int32 = 10
    Nn::Int32 = 10
    Nrp::Int32 = 10
    Nrn::Int32 = 10
    Model::String  = "SPM"
    Time::Array{Float64} = [0 3600]
    MeshType::String  = "L2"
    gsOrder::Int32 = 2
    Dimension::Int32 = 1
    #opt.Load = {"constant discharge 1C for 1h"}
    Current::Function = x-> 0
    CoupleMethod:: String  = "fully coupled"
    CoupleOrder::Int32 = 0
    y0::Array{Float64} = []
    dt::Array{Float64} = [0.1, 0.1]
    dtType::String  = "constant" # auto or manual
    dtThreshold::Float64 = 0.01
    SolveType::String  = "Crank-Nicolson" # forward, backward or Crank-Nicolson
    OutputType::String  = "auto" # auto or manual
    OutputTime::Array{Float64,1} = []
    Jacobi::String = "constant" # constant or update
end