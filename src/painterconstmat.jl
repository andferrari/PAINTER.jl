###################################################################################
# Antony Schutz 2015, ANR - POLCA
###################################################################################
# Pre calculus, MATRIX Inversion - Method involving orthogonal matrices
# ---------------------------------------------------------------------------------
function invmat_par(F3D::Array,rho_y::Real,eta::Real,nw::Int)#,paral::Bool)
# step IV of [0]
# inverse of inner matrix of C^\lambda_n
# [F F^H + eta/rho I]^-1
# Compute the inverse matrix of size nb*nb for image reconstruction, pre calculus step
# F3D: non uniform 3D Fourier Transform Matrix
# rho_x: admm parameters
# eta: (nwvlt*OIDATA.rho_spat+OIDATA.rho_spec+OIDATA.rho_ps+OIDATA.epsilon)
# if paral
mat2inv  = {(F3D[n],(rho_y/eta)) for n =1:nw }
M        = pmap(toinv,mat2inv)
# else
# # # # SERIAL
# nb = size(F3D[1],1)
# UnNb        = speye(nb)
# Mtmp = zeros(Complex{Float64},nb,nb,nw)
# M = { Mtmp[:,:,n] for n =1:nw }
# for n       = 1 : nw
#   M[n] = inv(rho_y/eta*F3D[n]*F3D[n]' + UnNb)
# end
# end
  return M
end
# ---------------------------------------------------------------------------------
# for parallel calculus of inverse
function toinv(mat2inv)
inv(mat2inv[2]*mat2inv[1]*mat2inv[1]' + speye(size(mat2inv[1],1)))
end

# ---------------------------------------------------------------------------------
# Phases To Phases difference Matrix
# ---------------------------------------------------------------------------------
# Method 1
function phasetophasediff(Closure_index::Matrix,nw::Integer,nb::Integer,T3::Integer,DP::Integer)
# create the matrix which link the phases to the phases difference
# Eq. 20 of [1] , 8 of [0]
# Closure phase and differential phase to phase matrix
# closure_index: index of phase closure (index related to phase vector index)
# nb is the number of base
# nw is the number of wavelength
# T3 and DP TBD for later, inform if TP and DP are present in data
# defined as in equation 14 of PAINTER
  nbnw    = nb*nw
if T3 == 1
  NT3     = size(Closure_index,1)
  rowt    = repeat([1:NT3*nw],inner=[3])
  colt    = repeat(vec(Closure_index'),outer=[nw]) + repeat([0:nw-1]*nb,inner=[3*NT3])
  Valuet  = repeat([1,1,-1],outer=[NT3*nw])
end
# Differential phase to phase matrix
# defined as in equation 19 of PAINTER [1], 7 of [0]
# rank = Nbas*(Nwvl-1)
if DP ==1
  rowb    = repeat([1:nb*(nw-1)],inner=[2])
  colb    = vec([repmat([1:nb]',1,nw-1),[nb+1:nbnw]'])
  Valueb  = repeat([1,-1],outer=[nb*(nw-1)])
end
# phase difference to phase matrix
# defined as in equation 20 of PAINTER [1], 6 of [0]
if(T3==1&&DP==0)
  HT3     = sparse(rowt,colt,Valuet,NT3*nw,nbnw)
  return HT3
elseif(T3==0&&DP==1)
  HDP     = sparse(rowb,colb,Valueb,nb*(nw-1),nbnw)
  return HDP
elseif(T3==1&&DP==1)
  HT3     = sparse(rowt,colt,Valuet,NT3*nw,nbnw)
  HDP     = sparse(rowb,colb,Valueb,nb*(nw-1),nbnw)
  H       = vcat(HT3,HDP)
  return H
end
end
###################################################################################
# Tools for NUFFT
# ---------------------------------------------------------------------------------
function planarray_par(tab_u::Array,tab_v::Array,nx::Int,nw::Int)#,parral::Bool)
# create Array for plan for non uniform fft
# tab_u and tab_v are spatial frequencies
# nx is images side size
# nw is the number of wavelength
# if parral
Mtmp = {((hcat(tab_u[:,n]/nx,tab_v[:,n]/nx))',nx) for n = 1 : nw}
return pmap(toplan,Mtmp)
# else
# # # SERIAL
# for n       = 1 : size(tab_u,2)
# Freq        = (hcat(tab_u[:,n]/nx,tab_v[:,n]/nx))'
# plan[n]     = NFFTPlan(Freq,(nx,nx))
# end
# end
end
# ---------------------------------------------------------------------------------
# for parallel calculus
function toplan(Mtmp)
NFFTPlan(Mtmp[1],(Mtmp[2],Mtmp[2]))
end
###################################################################################
# Matrices involved in data model
# ---------------------------------------------------------------------------------
# Non Uniform Discrete Fourier Matrices for Polychromatic data
# ---------------------------------------------------------------------------------
function nudft3d_par(tab_u::Array,tab_v::Array,nx::Int,nw::Int)
# create a list of non uniform dft Matrix per channel
# 3D nudft [number of bases (spatiale frequencies) per wavelength * number of pixels * number of wavelength]
# tab_u and tab_v must be of the form : [number of bases * number of wavelength]
# nx is images side size
# nw is the number of wavelength
# defined as in equation 3 of PAINTER [1]
if nx<1||nx<1
  error("nudft3d:Image size must be positive")
end
  nfu  	= length(tab_u)
  nfv  	= length(tab_v)
  if(nfu==nfv&&nfu>0&&nfv>0)
    # test if U and V have same size and have same number of wavelength, otherwise PB
    size_tab_u = size(tab_u)
    size_tab_v = size(tab_v)
    if(size_tab_u[2]==size_tab_v[2])
      (nb,nw)= size_tab_u
	  # # # SERIAL
	  #       F3D = complex(zeros(nb,nx*nx,nw))
	  # tic()
	  #       for n      = 1:nw
	  #         nudftmat   = non_uniform_dft(tab_u[:,n],tab_v[:,n])
	  #         F3D[:,:,n] = nudftmat
	  #       end
	  # toc()
	  # tic()
	  UVMAT  = {(tab_u[:,n],tab_v[:,n],nx,nb) for n =1:nw }
	  F3D    = pmap(non_uniform_dft_par,UVMAT)
	  # toc()
    return F3D
    else
    error("nudft3d: U ad V vectors have not of same number of wavelength")
    end
  else
    error("nudft3d: U ad V vectors are not of same size")
  end
end
# ---------------------------------------------------------------------------------
function non_uniform_dft_par(UVMAT)
# Non Uniform DFT Matrix:
# UVMAT [u,v,nx]
# for a nx*nx pixels image
# u, v vector of normalized spatial frequencies (between -0.5 and 0.5)
u 		= vec(UVMAT[1])
v 		= vec(UVMAT[2])
nx 		= UVMAT[3]
nb 		= UVMAT[4]

kx 		= [-nx/2:-1+nx/2]
k1  	= repmat(kx,1,nx)
k2  	= k1'
k1v 	= vec(k1)'/nx
k2v 	= vec(k2)'/nx
k1m 	= repmat(k1v,nb,1)
k2m 	= repmat(k2v,nb,1)

um 		= repmat(u,1,nx*nx)
vm 		= repmat(v,1,nx*nx)

return exp(-2im*pi*(k1m.*um+k2m.*vm))./nx
end
