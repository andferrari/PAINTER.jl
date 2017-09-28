###################################################################################
# Antony Schutz 2015, ANR - POLCA - 2016
###################################################################################
# Pre calculus, MATRIX Inversion - Method involving orthogonal matrices
# ---------------------------------------------------------------------------------
function invmat_par(F3D::Array,rho_y::Real,eta::Real,nw::Int)
# step IV of [0]
# inverse of inner matrix of C^\lambda_n
# [F F^H + eta/rho I]^-1
# Compute the inverse matrix of size nb*nb for image reconstruction, pre calculus step
# F3D: non uniform 3D Fourier Transform Matrix
# rho_x: admm parameters
# eta: (nwvlt*OIDATA.rho_spat+OIDATA.rho_spec+OIDATA.rho_ps+OIDATA.epsilon)
    mat2inv = [(F3D[n], (rho_y / eta)) for n in 1:nw ]
    M = pmap(toinv, mat2inv)
    return M
end
# ---------------------------------------------------------------------------------
# for parallel calculus of inverse
function toinv(mat2inv)
   return inv((mat2inv[2] * mat2inv[1] * mat2inv[1]') + speye(size(mat2inv[1], 1)))
end
# ---------------------------------------------------------------------------------
# Von Mises Weight estimation
# ---------------------------------------------------------------------------------
function ItKappa(K::Array{Float64,1};pre=1e-9)
   costK(Km,kn) = abs2( Km - (1 - besseli(1,kn)./besseli(0,kn) ) )
    NK = length(K)
    est = zeros(NK)
    varK = K ./ pi
    for m in 1:NK
        k,i = IterInKappa(costK,varK[m],0.,500.)
        it = 0
        while (sqrt(costK(varK[m],k[i]))>pre)&&(it<10)
            it+=1
            k,i = IterInKappa(costK,varK[m],k[max(i-1,1)],k[min(i+1,length(k))],N=100)
        end
        est[m] = k[i]
    end
    return est
end
function IterInKappa(cost::Function, K::Float64,a::Float64,b::Float64;N=1000)

    Kappa = collect(linspace(a,b,N))
    est = 0.
    tmp = zeros(N)
    for n in 1 : N
        tmp[n] = cost(K,Kappa[n])
    end
    return Kappa, indmin(tmp)
end

# ---------------------------------------------------------------------------------
# Phases To Phases difference Matrix
# ---------------------------------------------------------------------------------
# Method 1
function phasetophasediff(Closure_index::Matrix,nw::Integer,nb::Integer,T3::Integer,DP::Integer,dptype::String,dpprm::Integer)
# create the matrix which link the phases to the phases difference
# Eq. 20 of [1] , 8 of [0]
# Closure phase and differential phase to phase matrix
# closure_index: index of phase closure (index related to phase vector index)
# nb is the number of base
# nw is the number of wavelength
# T3 and DP TBD for later, inform if TP and DP are present in data
# defined as in equation 14 of PAINTER
# dptype is the kind of differential phases: all, ref, diag, frame, sliding
# dpprm  is the parameters associated to ref (lambda ref) or frame and sliding (horizon)

    nbnw = (nb * nw)
    if(T3 == 1)
      NT3 = size(Closure_index, 1)
      rowt = repeat(collect(1:(NT3 * nw)), inner = [3])
      colt = repeat(vec(Closure_index'), outer = [nw]) + repeat(collect(0:(nw - 1)) * nb, inner = [3 * NT3])
      Valuet = repeat([1, 1, -1], outer=[NT3 * nw])
    end

# Differential phase to phase matrix
# defined as in equation 19 of PAINTER [1], 7 of [0]
# rank = Nbas*(Nwvl-1)
    if(DP == 1)
      if dptype == "all"
          # rank = Nbas*( Nwvl - 1 )
          HDP = alldp(nb, nw)
      elseif dptype == "ref"
          # rank = Nbas*( Nwvl - 1 )
          # dpprm is the reference channel
          HDP = alldpref(nb, nw, dpprm)
      elseif dptype == "diag"
          # rank = Nbas*( Nwvl - 1 )
          HDP = diagdp(nb, nw)
      elseif dptype == "frame"
          # rank = Nbas*floor( Nwvl - Nwvl/horizon )
          # dpprm is the horizon size of window
          HDP = framedp(nb, nw, dpprm)
      elseif dptype == "sliding"
          # rank = Nbas*( Nwvl - 1 )
          # dpprm is the horizon size of window
          HDP = slidingdp(nb, nw, dpprm)
      elseif dptype == "phase"
          # rank = Nbas*Nwvl (full)
          HDP = nodp(nb, nw)
      else
          error("wrong choice for Differential Phase to Phase Matrix: all (default), frame, ref, diag, sliding or 'phase' ")
      end
    end

# phase difference to phase matrix
# defined as in equation 20 of PAINTER [1], 6 of [0]
    if((T3 == 1)&&(DP == 0))
      HT3 = sparse(rowt, colt, Valuet, (NT3 * nw), nbnw)
      return HT3

    elseif((T3 == 0)&&(DP == 1))
      return HDP

    elseif((T3 == 1)&&(DP == 1))
      HT3 = sparse(rowt, colt, Valuet, NT3 * nw, nbnw)
      H = vcat(HT3, HDP)
      return H
    end
end
# All differential phase, Painter V1
function alldp(nb::Int64,nw::Int64)
    nbnw = (nb * nw)
    rowb = repeat(collect(1:(nb * (nw - 1))), inner = [2])
    colb = vec(vcat(repmat(collect(1:nb)', 1, (nw-1)),collect((nb + 1):nbnw)'))
    Valueb = repeat(vcat(1,-1), outer=[nb*(nw-1)])
    HDP = sparse(rowb, colb, Valueb, nb * (nw-1), nbnw)
    return HDP
end
# case visphi is phases of complexe visibilities
function nodp(nb::Int64,nw::Int64)
    nbnw = (nb * nw)
    HDP = speye(nbnw)
    return HDP
end
# All differential phase, with reference channel defined by lambdaref
function alldpref(nb::Int64,nw::Int64,lambdaref::Int64)
    tmp = (lambdaref-1)*nb + collect(1:nb)'
    nbnw = (nb * nw)
    rowb = repeat(collect(1:(nb * (nw - 1))), inner = [2])

    colb = vec( vcat( repmat(tmp, 1, nw - 1), hcat( collect( 1:(tmp[1] - 1) )', collect( (1 + tmp[end]):nbnw )' ) ))
    Valueb = repeat(vcat(1,-1), outer=[nb*(nw-1)])
    HDP = sparse(rowb, colb, Valueb, nb * (nw-1), nbnw)
   return HDP
end
# All differential phase, Painter V1
function diagdp(nb::Int64,nw::Int64)
    nbnw = (nb * nw)
    tmp = collect(1:(nb * (nw - 1)))
    rowb = vcat(tmp,tmp)
    colb = vcat(tmp,collect((nb + 1):nbnw))
    Valueb = repeat(vcat(1,-1), inner=[nb*(nw-1)])
    HDP = sparse(rowb, colb, Valueb, nb * (nw-1), nbnw)
    return HDP
end

# Non overlapping windows for DP, Painter SPIE 2016
function framedp(nb::Int64,nw::Int64,horiz::Int64)
    horiz = horiz-1;
    nbloc   = floor(nw/(horiz+1));
    rowb    = Int64[]
    Colb    = Int64[]
    Valueb  = Int64[]
    for n   in 1:nbloc
        rowbtm = (n-1)*nb*(horiz) + repeat(collect(1:(nb * horiz)), inner = [2])
        Colbtm = (n-1)*nb*(horiz+1) + vcat( repmat(collect(1:nb)',1,horiz), collect(nb+1:nb*(horiz+1))' )
        Valuebtm = repmat(vcat(1,-1),nb*horiz,1)
        rowb    = vcat(rowb, round(Int,rowbtm))
        Colb    = vcat(Colb, vec(round(Int,Colbtm)))
        Valueb  = vcat(Valueb, vec(round(Int,Valuebtm)))
    end
    # % rest
    horiz2  = round(Int, (nb*nw-maximum( vec(Colb) ))/nb-1 )
    if horiz2>0
        rowbtm  = round(Int, maximum( vec(rowb) )) +repeat(collect(1:(nb * horiz2)), inner = [2])
        Colbtm  = round(Int, maximum( vec(Colb) )) + vcat( repmat(collect(1:nb)',1,horiz2), collect(nb+1:nb*(horiz2+1))' )
        Valuebtm = repmat(vcat(1,-1),nb*horiz2,1)
        rowb    = vcat(rowb, rowbtm)
        Colb    = vcat(Colb, vec(Colbtm))
        Valueb  = vcat(Valueb, vec(Valuebtm))
    end
    H       = sparse(rowb, Colb, Valueb, maximum(rowb),nb*nw)
    return H
end

# sliding windows, overlap 1, for DP, Painter SPIE 2016
function slidingdp(nb::Int64,nw::Int64,horiz::Int64)
    horiz = horiz-1;
    nbnw = (nb * nw)
    row1 = repeat(collect(1:(nb*horiz)), inner = [2])
    col1 = vec(vcat(repmat(collect(1:nb)', 1, horiz),collect((nb + 1):(nb*(horiz+1)))'))
    value1 = repeat(vcat(1,-1), outer=[nb*horiz])

    row2    = Int64[]
    col2    = Int64[]
    value2  = Int64[]

    vnb = collect(1:nb)
    for n   in 1:(nw-1-horiz)
        rowbtm = (n-1)*nb + maximum(row1) + repeat(vnb, inner = [2])
        colbtm = vcat(n*nb + vnb' , (n-1)*nb + maximum(col1) + vnb')
        valuebtm = repmat(vcat(1,-1),nb,1)
        row2 = vcat(row2, round(Int,rowbtm))
        col2 = vcat(col2, vec(round(Int,colbtm)))
        value2 = vcat(value2, vec(round(Int,valuebtm)))
    end

    row = vcat(row1, row2)
    col = vcat(col1, vec(col2))
    value = vcat(value1, vec(value2))
    H = sparse(row, col, value, nb*(nw-1), nb*nw)
    return H
end

###################################################################################
# Tools for NUFFT
# ---------------------------------------------------------------------------------
function planarray_par(tab_u::Array,tab_v::Array,nx::Int,nw::Int)
# create Array for plan for non uniform fft
# tab_u and tab_v are spatial frequencies
# nx is images side size
# nw is the number of wavelength
    Mtmp = [( hcat(tab_u[:,n] / nx, tab_v[:,n] / nx)', nx) for n in 1:nw]
    return pmap(toplan, Mtmp)
end
# ---------------------------------------------------------------------------------
# for parallel calculus
function toplan(Mtmp)
    NFFTPlan(Mtmp[1], (Mtmp[2], Mtmp[2]))
end
###################################################################################
# Matrices involved in data model
# ---------------------------------------------------------------------------------
# Non Uniform Discrete Fourier Matrices for Polychromatic data
# ---------------------------------------------------------------------------------
function nudft3d_par(tab_u::Array,tab_v::Array,nb::Int,nx::Int,nw::Int)
# create a list of non uniform dft Matrix per channel
# 3D nudft [number of bases (spatiale frequencies) per wavelength * number of pixels * number of wavelength]
# tab_u and tab_v must be of the form : [number of bases * number of wavelength]
# nx is images side size
# nw is the number of wavelength
# defined as in equation 3 of PAINTER [1]
    if(nx < 1)
        error("nudft3d:Image size must be positive")
    end
    nfu = length(tab_u)
    nfv = length(tab_v)

    if((nfu == nfv)&&(nfu > 0)&&(nfv > 0))

# test if U and V have same size and have same number of wavelength, otherwise PB
        if(size(tab_u,2) == size(tab_v,2))
            UVMAT = [(tab_u[:, n], tab_v[:, n], nx, nb) for n in 1:nw ]
	          F3D = pmap(non_uniform_dft_par, UVMAT)
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
# UVMAT[1] : U
# UVMAT[2] : V
# UVMAT[3] : nx
# UVMAT[4] : nb
    k1  	= repmat(collect((-UVMAT[3] / 2):(-1 + (UVMAT[3] / 2))), 1, UVMAT[3])
    k1m 	= repmat((vec(k1 )' / UVMAT[3]), UVMAT[4], 1)
    k2m 	= repmat((vec(k1')' / UVMAT[3]), UVMAT[4], 1)
    um 		= repmat(vec(UVMAT[1]), 1, (UVMAT[3] * UVMAT[3]))
    vm 		= repmat(vec(UVMAT[2]), 1, (UVMAT[3] * UVMAT[3]))
    return exp.(-2 * im * pi * ((k1m .* um) + (k2m .* vm) )) ./ UVMAT[3]
end
