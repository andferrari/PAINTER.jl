###################################################################################
# Antony Schutz 2015, ANR - POLCA
###################################################################################
## Data Type
###################################################################################
# PAINTER Data Type
#
# Structure containing all user data and oifits data
# OIDATA::PAINTER_Input
#
type PAINTER_Input
  PlotFct::Function
  Folder::ASCIIString
  FilesName::Array{ASCIIString}
  indfile::Array{Int,1}
  indwvl::Array{Int,1}
  wvl::Array{Float64}
  U::Array{Float64}
  V::Array{Float64}
  P::Array{Float64}
  W::Array{Float64}
  Xi::Array{Float64}
  K::Array{Float64}
  Closure_index::Array{Int}
  nb::Int
  nw::Int
  nx::Int
  FOV::Real
  lambda_spat::Real
  lambda_spec::Real
  lambda_L1::Real
  rho_y::Real
  rho_spat::Real
  rho_spec::Real
  rho_ps::Real
  alpha::Real
  beta::Real
  eps1::Real
  eps2::Real
  epsilon::Real
  mask3D::Array
  xinit3D::Array
  Wvlt::Array
  paral::Bool
  T3::Array{Float64}
  T3err::Array{Float64}
  DP::Array{Float64}
  DPerr::Array{Float64}
end
# Structure containing all data which are modified during admm
# PDATA::PAINTER_Data
type PAINTER_Data
  eta::Real
  plan::Array
  F3D::Array
  H::SparseMatrixCSC
  M::Array
  #Wvlt::Array
  x::SharedArray{Float64}#Array
  vHt::Array
  z::Array
  Hx::Array
  tau_s::Array
  w::Array
  v::Array
  r::Array
  tau_w::Array
  tau_v::Array
  tau_r::Array
  Spcdct::Array
  Fx::SharedArray{Complex{Float64}}
  tau_xc::Array
  tau_pwc::Array
  tau_xic::Array
#   ys::Array
#   y_tampon::Array
  yc::Array
  y_v2::Array
  y_phi::Array
  crit1::Array
  crit2::Array
  ind::Int
  CountPlot::Int
  count::Int
  nbitermax::Int
end
# Structure for OptimPack, see OptimPack.jl
# https://github.com/emmt/OptimPack.jl
type OptOptions
  ls
  scl
  gat
  grt
  vt
  memsize
  mxvl
  mxtr
  stpmn
  stpmx
end

###################################################################################
# Initialise structure with nothing
###################################################################################
function optiminit(ls,scl,gat,grt,vt,memsize,mxvl,mxtr,stpmn,stpmx)
return OptOptions(ls,scl,gat,grt,vt,memsize,mxvl,mxtr,stpmn,stpmx)
end
function painterinputinit()
  return PAINTER_Input(painterplotfct,"",[],[],[],[],[],[],[],[],[],[],[],0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,[],[],[],true,[],[],[],[])
end
function painterdatainit()
  return PAINTER_Data(0.,[],[],speye(0),[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],Float64[],Float64[],0,0,0,0)
#   return PAINTER_Data(0.,[],[],speye(0),[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],Float64[],Float64[],0,0,0,0)
end
