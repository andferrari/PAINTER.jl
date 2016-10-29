###################################################################################
# Antony Schutz 2015, ANR - POLCA - 2016
###################################################################################
## Data Type
###################################################################################
# PAINTER Data Type
#
# Structure containing all user data and oifits data
# OIDATA::PAINTER_Input
#
###################################################################################
# Initialise structure with nothing
###################################################################################
type PAINTER_Input
    PlotFct::Function               # user defined plot function
    Folder::String             # Folder of OIFITS/FITS Files
    FilesName::Array{String}   # name of the files present in Folder
    indfile::Array{Int,1}           # index of file used
    indwvl::Array{Int,1}            # index of wavelength used
    wvl::Array{Float64}             # mean value of anayzed wavelength
    U::Array{Float64}               # U spatial frequencies
    V::Array{Float64}               # V spatial frequencies
    P::Array{Float64}               # Matrix of V2
    W::Array{Float64}               # Matrix of V2err
    Xi::Dict                        # Dictionary of Vector of phases difference
    K::Dict                         # Dictionary of Vector of phases difference error
    Closure_index::Array{Int}       # index of phases closure
    nb::Int                         # number of bases (spatiale frequencies per wavelength)
    nw::Int                         # number of wavelength
    nx::Int                         # side size of image in pixels
    FOV::Real                       # Field Of View in arcsecond
    lambda_spat::Real               # From here to end: admm variables
    lambda_spec::Real
    lambda_L1::Real
    rho_y::Real

    rho_y_gamma::Real
    rho_y_xi::Real

    rho_spat::Real
    rho_spec::Real
    rho_ps::Real
    alpha::Real
    beta::Real
    eps1::Real
    eps2::Real
    epsilon::Real
    mask3D::Array                   # support constraint array
    xinit3D::Array                  # initial estimate array
    Wvlt::Array                     # list of wavelets basis
    T3::Array{Float64}              # Matrix of phases closure
    T3err::Array{Float64}           # Matrix of phases closure error
    DP::Array{Float64}              # Matrix of Differential phases
    DPerr::Array{Float64}           # Matrix of Differential phases error
    isDP::Int                       # To know if there is DP in data
    dptype::String             # way to construct diff phases
    dpprm::Int                      # horizon for sliding dp or index of reference
    baseNb::Dict
    orderedCluster::Dict
    flux::Array
    norm::Real
    # Kbefore::Dict
end
# Structure containing all data which are modified during admm
# PDATA::PAINTER_Data
type PAINTER_Data
    eta::Real
    plan::Array
    F3D::Array
    H::Dict
    M::Array
    x::Array
    z::Array
    tau_s::Array
    w::Array
    v::Array
    r::Array
    tau_w::Array
    tau_v::Array
    tau_r::Array
    Fx::Array
    tau_xc::Array
    tau_pwc::Array
    tau_xic::Array
    yc::Array
    y_v2::Array
    y_phi::Array
    crit1::Array
    crit2::Array
    ind::Int
    CountPlot::Int
    count::Int

end

function painterinputinit()
    return PAINTER_Input( painterplotfct, "", [], [], [], [], [], [], []
                         , [], Dict{}(), Dict{}(), [], 0., 0., 0., 0., 0., 0., 0., 0.
                         , 0., 0.
                         , 0., 0., 0., 0., 0., 0., 0., 0., [], [], []
                         , [], [], [], [], 0, "", 0, Dict{}(), Dict{}(),[],0) # [],Dict{}())
end
function painterdatainit()

  return PAINTER_Data(0., [], [], Dict{}(), [], [], [], [], [], []
                      , [], [], [], [], [], [], [], [], [], []
                      , [], Float64[], Float64[], 0, 0, 1)
end
