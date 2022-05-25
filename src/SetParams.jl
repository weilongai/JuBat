"""
        Set up a structure for parameters of a lithium-ion cell
        The names of variables use the beginning uppercase words to specify the domain:    
            NE -- negative electrode
            PE -- positive electrode
            SP -- separator
            EL -- electrolyte
            NCC -- negative current collector
            PCC -- positive current collector
            Cell -- cell 
            Tab -- tab
            BD -- binder

        The variables are for properties of :
            theta -- stoichiometry
            thickness -- thickness [m]
            lambda -- thermal conductivities [W/(m K)]
            Ds -- solid phase diffusion coefficients [m^2/s]
            rho -- density [kg/m^3]
            heat_Q -- specific heat capacities [J/(kg K)]
            sig -- conductivity [S/m]
            eps -- porosity
            brugg -- Bruggeman coefficients
            k -- reaction rate constant [ A m^2.5/mol^1.5]
            tplus -- transference number 
            h -- heat exchange coefficient [W/(m^2 K)]
            Rs -- electrode particle radius [m]
            cs_max -- maximum concentration in solid phase [mol/m^3]
            ce0 -- electrolyte lithium-ions initial concentration [mol/m^3]
            U -- open-circuit potential (OCP) [V]
            dUdT -- entropic change of the OCP [V/K]
            theta_0 --  Theta @ 0% Lithium Concentration 
            theta_100 --  Theta @ 100% Lithium Concentration 
            alpha -- Alpha Factor
"""
@with_kw mutable struct Electrode
    theta_100::Float64 = 0
    theta_0::Float64 = 0
    thickness::Float64 = 0
    lambda::Float64 = 0
    Ds::Float64 = 0
    rho::Float64 = 0
    heat_Q::Float64 = 0
    eps::Float64 = 0
    eps_fi::Float64 = 0
    brugg::Float64 = 0
    k::Float64 = 0
    cs_max::Float64 = 0
    cs_0::Float64 = 0
    Rs::Float64 = 0
    as::Float64 = 0
    sig::Float64 = 0
    Eac_D::Float64 = 0
    Eac_k::Float64 = 0
    alpha::Float64 = 0
    OCP::Function = x-> 0.0
    dOCP::Function = x-> 0.0
end

@with_kw mutable struct Separator
    thickness::Float64 = 0
    lambda::Float64 = 0
    rho::Float64 = 0
    heat_Q::Float64 = 0
    eps::Float64 = 0
    eps_fi::Float64 = 0
    brugg::Float64 = 0
end

@with_kw mutable struct CurrentCollector
    thickness::Float64 = 0
    lambda::Float64 = 0
    rho::Float64 = 0 
    heat_Q::Float64 = 0
    sig::Float64 = 0 
end

@with_kw mutable struct Electrolyte
    De::Function = x-> 0
    rho::Float64 = 0 
    heat_Q::Float64 = 0
    tplus::Float64 = 0
    ce0::Float64 = 0
end

@with_kw mutable struct Cell
    length::Float64 = 0
    width::Float64 = 0
    # wrapper::Float64
    I1C::Float64 = 0
    no_layers::Int32 = 0 
    capacity::Float64 = 0 
    Cooling_surface::Float64 = 0 
    Total_surface :: Float64 = 0
    v_h::Float64 = 0 
    v_l::Float64 = 0
    volume::Float64 = 0 
    rho::Float64 = 0
    alphaT::Float64 = 0
    heat_Q::Float64 = 0 
    h::Float64 = 0
    T_ref::Float64 = 0
end

@with_kw mutable struct Tab
    length::Float64 = 0
    width::Float64 = 0
    area::Float64 = 0
end
# param_dim.Tab.width = 0.75 * param_dim.Tab.length;  
# param_dim.Tab.area = param_dim.Tab.length * param_dim.Tab.width;

@with_kw mutable struct Binder
    rho::Float64 = 0
end

@with_kw mutable struct Scale
    L::Float64 = 1e-6
    a0::Float64 = 1e6
    r0::Float64 = 1e-6
    t0::Float64 = 3600
    T0::Float64 = 298
    F::Float64 = 96485.33289
    R::Float64 = 8.314
    j::Float64 = 0
    Ds_p::Float64 = 0
    Ds_n::Float64 = 0
    ts_p::Float64 = 0
    ts_n::Float64 = 0
    te::Float64 = 0
    De::Float64 = 0
    phi::Float64 = 0
    sig::Float64 = 0
    kappa::Float64 = 0
    cp_max::Float64 = 0
    cn_max::Float64 = 0
    ce::Float64 = 0
    k_p::Float64 = 0
    k_n::Float64 = 0
    I_typ::Float64 = 0
end

@with_kw mutable struct Params
    PE::Electrode
    NE::Electrode
    EL::Electrolyte
    SP::Separator
    cell::Cell
    PCC::CurrentCollector
    NCC::CurrentCollector
    tab::Tab
    binder::Binder
    scale::Scale
end

function ChooseCell(CellType::String="LG M50")
"""
    This is a function choose a cell
    Input - CellType::String, including options of 
        1. "LG M50" for the LG M50 cells
    Output - param_dim::Params for all param_dim parameters
"""
    if CellType == "LG M50"
        include("../src/data/LGM50.jl") # pathof(JuBat)
    end
    param_dim.PE.as = 3 * (1 - param_dim.PE.eps) / param_dim.PE.Rs
    param_dim.NE.as = 3 * (1 - param_dim.NE.eps) / param_dim.NE.Rs
    param_dim.scale.I_typ = param_dim.cell.I1C / param_dim.cell.length / param_dim.cell.width / param_dim.cell.no_layers
    param_dim.scale.L = param_dim.PE.thickness + param_dim.NE.thickness + param_dim.SP.thickness;
    param_dim.scale.j = param_dim.scale.I_typ / param_dim.scale.a0 / param_dim.scale.L;
    param_dim.scale.Ds_p = param_dim.scale.r0 / param_dim.scale.F / param_dim.PE.cs_max * param_dim.scale.j;
    param_dim.scale.Ds_n = param_dim.scale.r0 / param_dim.scale.F / param_dim.NE.cs_max * param_dim.scale.j;
    param_dim.scale.ts_p = param_dim.scale.F * param_dim.PE.cs_max * param_dim.scale.r0 / param_dim.scale.j / param_dim.scale.t0;
    param_dim.scale.ts_n = param_dim.scale.F * param_dim.NE.cs_max * param_dim.scale.r0 / param_dim.scale.j / param_dim.scale.t0;
    param_dim.scale.te = param_dim.scale.F * param_dim.EL.ce0 / param_dim.scale.a0 / param_dim.scale.j / param_dim.scale.t0;
    param_dim.scale.De = param_dim.scale.a0 * param_dim.scale.L^2 / param_dim.scale.F / param_dim.EL.ce0 * param_dim.scale.j;
    param_dim.scale.phi = param_dim.scale.T0 * param_dim.scale.R / param_dim.scale.F;
    param_dim.scale.sig = param_dim.scale.L * param_dim.scale.I_typ / param_dim.scale.phi;
    param_dim.scale.kappa = param_dim.scale.L^2 * param_dim.scale.a0;
    param_dim.scale.cp_max = param_dim.PE.cs_max;
    param_dim.scale.cn_max = param_dim.NE.cs_max;
    param_dim.scale.ce = param_dim.EL.ce0;
    param_dim.scale.k_p = param_dim.scale.j / param_dim.PE.cs_max / sqrt(param_dim.EL.ce0);
    param_dim.scale.k_n = param_dim.scale.j / param_dim.NE.cs_max / sqrt(param_dim.EL.ce0);
    return param_dim
end

function NormaliseParam(param_dim::Params)
    """
        This is a function to normalise the parameters
        Input - param_dim::Params (with units) for a cell
        Output - param::Params (normalised)
    """
    # normalise the parameters, while thermal parameters are not covered, need revisit
    param = deepcopy(param_dim)

    # posotove electrode
    param.PE.theta_100 = param_dim.PE.theta_100
    param.PE.theta_0 = param_dim.PE.theta_0;
    param.PE.cs_0 = param_dim.PE.cs_0 / param_dim.PE.cs_max;
    param.PE.thickness = param_dim.PE.thickness / param.scale.L; 
    param.PE.Ds = param_dim.PE.Ds / param.scale.Ds_p;
    param.PE.eps = param_dim.PE.eps;
    param.PE.eps_fi = param_dim.PE.eps_fi;
    param.PE.brugg = param_dim.PE.brugg;
    param.PE.k = param_dim.PE.k / param.scale.k_p;
    param.PE.Rs = param_dim.PE.Rs / param.scale.r0;
    param.PE.sig = param_dim.PE.sig / param.scale.sig;
    param.PE.OCP = x-> param_dim.PE.OCP(x) / param.scale.phi;
    param.PE.as = param_dim.PE.as / param.scale.a0

    # negative electrode
    param.NE.theta_100 = param_dim.NE.theta_100;
    param.NE.theta_0 = param_dim.NE.theta_0;
    param.NE.cs_0 = param_dim.NE.cs_0 / param_dim.NE.cs_max;
    param.NE.thickness = param_dim.NE.thickness / param.scale.L;
    param.NE.Ds = param_dim.NE.Ds / param.scale.Ds_n;
    param.NE.eps = param_dim.NE.eps;
    param.NE.eps_fi = param_dim.NE.eps_fi;
    param.NE.brugg = param_dim.NE.brugg;
    param.NE.k = param_dim.NE.k / param.scale.k_n;
    param.NE.Rs = param_dim.NE.Rs / param.scale.r0;
    param.NE.sig = param_dim.NE.sig / param.scale.sig;
    param.NE.OCP =x-> param_dim.NE.OCP(x) / param.scale.phi;
    param.NE.as = param_dim.PE.as / param.scale.a0

    # separator
    param.SP.thickness = param_dim.SP.thickness / param.scale.L;
    param.SP.eps = param_dim.SP.eps;
    param.SP.eps_fi = param_dim.SP.eps_fi;
    param.SP.brugg = param_dim.SP.brugg;

    # positive current colloctor
    param.PCC.thickness = param_dim.PCC.thickness / param.scale.L;
    param.PCC.sig =  param_dim.PCC.sig / param.scale.sig;

    # negative current colloctor
    param.NCC.thickness = param_dim.NCC.thickness / param.scale.L;
    param.NCC.sig = param_dim.NCC.sig / param.scale.sig;

    # electrolyte
    param.EL.De = x-> param_dim.EL.De(x) / param.scale.De;
    param.EL.tplus = param_dim.EL.tplus;
    param.EL.ce0 = param_dim.EL.ce0 / param.scale.ce;

    return param
end

