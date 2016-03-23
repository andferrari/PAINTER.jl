###################################################################################
# Antony Schutz 2015, ANR - POLCA - 2016
###################################################################################
# REF
###################################################################################
# [0] PAINTER
# Schutz, A., Ferrari, A., Mary, D., Soulez, F., Thiébaut, E., Vannier, M.
# Large scale 3D image reconstruction in optical interferometry
# EUSIPCO 2015
#
# [1] PAINTER
# Schutz, A., Ferrari, A., Mary, D., Soulez, F., Thiébaut, E., Vannier, M.
# PAINTER: a spatiospectral image reconstruction algorithm for optical interferometry
# J. Opt. Soc. Am. A, Vol. 31, N. 11, pp 2334--2345, 2014.
#
# [2] ADMM
# S. Boyd, N. Parikh, E. Chu, B. Peleato and J. Eckstein.
# “Distributed optimization and statistical learning
# via the alternating direction method of multipliers,”
# Found. Trends Mach. Learn., 3(1):1–122, January 2011.
###################################################################################
# # Operation involved in PAINTER
# ---------------------------------------------------------------------------------
# # # # # # # # # # #            W A V E L E T S
# ---------------------------------------------------------------------------------
# Object Estimation - Orthognal matrix - wavelet
# ---------------------------------------------------------------------------------
function estimx_par{Tw<:WT.OrthoWaveletClass}(x::SharedArray{Float64,3},Fx::SharedArray{Complex{Float64},2},
    rho_y::Float64,rho_spat::Float64,rho_spec::Float64,rho_ps::Float64,eta::Float64,
    yc::Array{Complex{Float64},2},z::Array{Float64,4},v::Array{Float64,3},w::Array{Float64,3},
    tau_xc::Array{Complex{Float64},2},tau_s::Array{Float64,4},tau_v::Array{Float64,3},tau_w::Array{Float64,3},
    nb::Int,nw::Int,nx::Int,NWvlt::Int,plan::Array{Any,1},Wvlt::Array{Tw,1},M::Array{Any,1},paral::Bool)
# Estimate the constrained, regularized 3D images from complexe visibilities
# step IV of PAINTER [0]
#
# x: images [(nx*nx)*nw], Fx is nb*nw is Non uniform Fourier transform of 3D images
# rho_x and tau_x: admm parameters and Lagrange mutlipliers
# yc: estimated complexe visibilities
# plan is the non uniform fft plan
# Wvlt is the list of used wavelets basis
# M is the used inverse matrix pre computed
    maty = (rho_y * yc) - tau_xc
    matz = (rho_spat * z) - tau_s
    matv = (rho_spec * v) - tau_v
    matw = (rho_ps * w) - tau_w
    Reg  = matv + matw

# parallel sum of wavelet basis
    if paral
        MAP = [(matz[:, :, n, b] , wavelet(Wvlt[b]) ) for n in 1:nw, b in 1:NWvlt]
        wvd = sum(reshape(pmap(myidwt,MAP), nw, NWvlt ), 2)

# parallel image reconstruction
        @sync @parallel for n in 1:nw
            xtmp = (nfft_adjoint(plan[n], maty[:, n]) / nx) + Reg[:, :, n] + wvd[n]
            xfst = nfft_adjoint(plan[n], M[n] * (nfft(plan[n], xtmp) / nx)) / nx
            xtmp = (xtmp - (rho_y / eta) * xfst) / eta
            x[:,:,n] = real(xtmp)
            Fx[:,n]  = nfft(plan[n], xtmp) / nx
        end

    else
# # SERIAL
        wavdec = zeros(nx, nx, nw, nb)
        for b in 1:NWvlt, n in 1:nw
            wavdec[:, :, n, b] = idwt(matz[:, :, n, b], wavelet(Wvlt[b]))
        end
        wvd = sum(wavdec, 4)
        for n in 1:nw
            xtmp = (nfft_adjoint(plan[n], maty[:, n]) / nx) + Reg[:, :, n]  + wvd[:, :, n]
            xfst = nfft_adjoint(plan[n], M[n] * (nfft(plan[n], xtmp) / nx)) / nx
            xtmp = (xtmp - (rho_y / eta) * xfst) / eta
            x[:,:,n] = real(xtmp)
            Fx[:,n]  = nfft(plan[n], xtmp) / nx
        end
    end
    return x,Fx
end
# ---------------------------------------------------------------------------------
# for parallel calculus of inverse wavelet basis
function myidwt(M)
   return idwt(M[1], M[2])
end
###################################################################################
# Regularization - constraint
# Section 4.B of PAINTER
###################################################################################
# Proximal operator for V2 and Cardano's formula
# ---------------------------------------------------------------------------------
# estimate the complexe visibilities from the squared modulus visibilities
# proximal operator
# P matrix of observed V2 [Nb*Nwvl]
# W matrix of observed V2 error variance [Nb*Nwvl]
# rho_y: Admm parameter (scalar>0)
# alpha: relative weight compared to phases difference: (scalar>0)
# (Global Cost) = alpha * (cost of V2) + Beta * (Cost of Phases diffrence)
# y_v2: auxiliary variable (last estimate of complexe visibility)
# nb is the number of base
# nw is the number of wavelength
#
# Equation 58-59 of PAINTER
#
function proxv2!(y_v2::Array,P::Array,rho_y::Real,alpha::Real,nb::Int,nw::Int)
    mod_y = abs(y_v2)
    ang_y = angle(y_v2)
    tmp1 = rho_y / (4 * alpha)
    tmp2 = alpha
    for m in 1:nb,n in 1:nw
# # without weight
        sol = max(0., cubicroots(1., 0., tmp1 -P[ind], -tmp1 * mod_y[m,n]))
        cst = tmp2*(P[ind] - sol.^2).^2 + .5 .* rho_y * (sol - mod_y[m,n] ).^2
        (a,b) = findmin(cst)
        mod_y[m,n] = sol[b]
    end
    y_v2[:] = mod_y .* exp(im .* ang_y)
end
function proxv2!(y_v2::Array,P::Array,W::Array,rho_y::Real,alpha::Real,nb::Int,nw::Int)
    mod_y = abs(y_v2)
    ang_y = angle(y_v2)
    tmp1 = W.*rho_y /(4 * alpha)
    tmp2 = alpha ./ W
    for m in 1:nb, n in 1:nw
        sol = max(0., cubicroots(1., 0., tmp1[m,n] - P[m,n], -tmp1[m,n] .* mod_y[m,n]))
        cst = tmp2[m,n] .* (P[m,n] - sol.^2).^2 + .5 .* rho_y * (sol - mod_y[m,n] ).^2
        (a,b) = findmin(cst)
        mod_y[m,n] = sol[b]
    end
    y_v2[:] = mod_y .* exp(im .* ang_y)
end
# ---------------------------------------------------------------------------------
# ----- Cardano's formula
function cubicroots(a::Real,b::Real,c::Real,d::Real)
# finds real valued roots of cubic
# Arguments:
#     a, b, c, d - coeffecients of cubic defined as
#                  ax^3 + bx^2 + cx + d = 0
# Returns:
# root   - an array of 1 or 3 real valued roots
# Reference:  mathworld.wolfram.com/CubicFormula.html
# Code follows Cardano's formula
# Copyright (c) 2008 Peter Kovesi
# School of Computer Science & Software Engineering
# The University of Western Australia
# pk at csse uwa edu au
# http://www.csse.uwa.edu.au/
# realcuberoot - computes real-valued cube root
function realcuberoot(x::Real)
    sign(x) .* abs(x).^(1 / 3)
end
    # Divide through by a to simplify things
    b = b / a
    c = c / a
    d = d / a
    bOn3  = b/3.
    q = (3 .* c - b^2) / 9.
    r = (9 .* b * c - 27 .* d - 2 * b^3) / 54.
    discriminant = q^3 + r^2
    if discriminant >= 0        # We have 1 real root and 2 imaginary
        s = realcuberoot(r + sqrt(discriminant))
        t = realcuberoot(r - sqrt(discriminant))
        root = s + t - bOn3     # Just calculate the real root
    else                        # We have 3 real roots
        # In this case (r + sqrt(discriminate)) is complex so the following
        # code constructs the cube root of this complex quantity
        rho = sqrt(r^2 - discriminant)
        cubeRootrho = realcuberoot(rho)    # Cube root of complex magnitude
        thetaOn3 = acos(r / rho) / 3       # Complex angle/3
        crRhoCosThetaOn3 = cubeRootrho * cos(thetaOn3)
        crRhoSinThetaOn3 = cubeRootrho * sin(thetaOn3)
        root = zeros(3)
        root[1] = 2 * crRhoCosThetaOn3 - bOn3
        root[2] = -crRhoCosThetaOn3 - bOn3 - sqrt(3) * crRhoSinThetaOn3
        root[3] = -crRhoCosThetaOn3 - bOn3 + sqrt(3) * crRhoSinThetaOn3
    end
  return root
end
# ---------------------------------------------------------------------------------
# Proximal operator for phases difference (Section 5.2 PAINTER)
# ---------------------------------------------------------------------------------
function proxphase!(y_phi::Matrix,Xi::Vector,K::Vector,rho_y::Real,beta::Real,nb::Int,nw::Int,H::SparseMatrixCSC,OPTOPT::OptOptions)
# estimate the complexe visibilities from the phases differences
# y_phi: vector of auxiliary variable \tilde{y} (last estimate of complexe visibility) is upadted
# Xi vector of observed Phases Difference [ {Nwvl*(Nb-1)*(Nb-2)/2 + Nb*(Nwvl-1)} *{1}]
# K the variance of phases must be corrected using Von Mises Ditribution (function EstKapVM)
# rho_y: Admm parameter (scalar>0)
# beta: relative weight compared to phases difference: (scalar>0)
# (Global Cost) = alpha * (cost of V2) + Beta * (Cost of Phases diffrence)
# nb is the number of base
# nw is the number of wavelength
# H: Phases difference to Phases matrix (function Ph2PhDiff)
# OPTOPT: structure of OptimPack vmlm option, see optimpack
    y_t = vec(y_phi)
    gam_t = abs(y_t)
    phi_t = angle(y_t)
    phi_0 = angle(y_t)
    function cost!{T<:Real}(x_phi::Array{T,1}, g_phi::Array{T,1})
        return costgradphi!(x_phi, g_phi, gam_t, phi_t, y_t, Xi, K, beta, rho_y, H)
    end

    phi = OptimPack.vmlm(cost!, phi_0, OPTOPT.memsize, verb = OPTOPT.vt, scaling = OPTOPT.scl
                         , grtol = OPTOPT.grt, gatol=OPTOPT.gat, lnsrch=OPTOPT.ls, maxeval=OPTOPT.mxvl
                         , maxiter=OPTOPT.mxtr, stpmin=OPTOPT.stpmn, stpmax=OPTOPT.stpmx)

    if phi!=nothing
        Ek = phi_t - phi
        gam = max(0.0, gam_t .* cos(Ek))
        y_phi[:] = reshape(gam .* exp(im * phi), nb, nw)
    else
        y_phi[:] = reshape(gam_t .* exp(im * phi_t), nb, nw)
    end
end
function costgradphi!(x_phi::Vector,g_phi::Vector,gam_t::Vector,phi_t::Vector,y_t::Vector,Xi::Vector,K::Vector,beta::Real,rho_y::Real,H::SparseMatrixCSC)
# cost to minimize for phase difference minimization
# x_phi: Phases vector to estimate
# g_phi: gradiant of phase estimation cost function
# gam_t modulus of y_tilde
# phi_t phases of y_tilde
# Xi vector of observed Phases Difference [ {Nwvl*(Nb-1)*(Nb-2)/2 + Nb*(Nwvl-1)} *{1}]
# K the Von mises variance
# rho_y: Admm parameter (scalar>0)
# beta: relative weight compared to phases difference: (scalar>0)
# (Global Cost) = alpha * (cost of V2) + Beta * (Cost of Phases diffrence)
# H: Phases difference to Phases matrix (function Ph2PhDiff)
    Ek = phi_t - x_phi
    gam = max(0.0, gam_t .* cos(Ek))
    dphi = H * x_phi - Xi;
    yest = gam .* exp(im * x_phi)
    w1 = -sum(K .* cos(dphi))
    w2 = sum(abs(yest - y_t).^2)
    g_phi[:]= beta .* H' * (sin(dphi) ./ K) - rho_y .* gam .* gam_t .* sin(Ek)
    f = ((beta * w1) + (rho_y * w2)) / 2
    return f
end
###################################################################################
# MAIN ADMM LOOP
###################################################################################
# method 1, the algorithm is already initialized, structures are created and filed
# the admm can work with this informations
# function painteradmm(PDATA::PAINTER_Data,OIDATA::PAINTER_Input,OPTOPT::OptOptions,nbitermax::Int,aff::Bool)
#     painteradmm(PDATA,OIDATA,OPTOPT,nbitermax,aff)
# end
# for information about parameters read PAINTER [1] and admm [2]
function painteradmm(PDATA::PAINTER_Data,OIDATA::PAINTER_Input,OPTOPT::OptOptions,nbitermax::Int,aff::Bool)
    const nx = OIDATA.nx
    const nb = OIDATA.nb
    const nw = OIDATA.nw
    const wvl = OIDATA.wvl
    const lambda_spat = OIDATA.lambda_spat
    const lambda_spec = OIDATA.lambda_spec
    const lambda_L1 = OIDATA.lambda_L1
    const epsilon = OIDATA.epsilon
    const rho_y = OIDATA.rho_y
    const rho_spat = OIDATA.rho_spat
    const rho_spec = OIDATA.rho_spec
    const rho_ps = OIDATA.rho_ps
    const alpha = OIDATA.alpha
    const beta = OIDATA.beta
    const eps1 = OIDATA.eps1
    const eps2 = OIDATA.eps2
    const Wvlt = OIDATA.Wvlt
    const mask3D = OIDATA.mask3D
    const P = OIDATA.P
    const W = OIDATA.W
    const Xi = OIDATA.Xi
    const K = OIDATA.K
    const paral = OIDATA.paral
    const eta = PDATA.eta
    const plan = PDATA.plan
    const F3D = PDATA.F3D
    const H = PDATA.H
    const M = PDATA.M
    const NWvlt = length(Wvlt)
# ----------------------------------
# Check if Pyplot is used to graphics
    if aff
# check if Pyplot is installed
        if(Pkg.installed("PyPlot") == nothing)
            println("")
            println("PyPlot is not installed: Pkg.add(''PyPlot''), aff=false ")
            aff = false
        else # PyPlot is installed, check if Painter know where is PyPlot
            try
                PyPlot.pygui(true)
            catch e
                println("")
                println("PyPlot not used: using PyPlot, try to load it otherwise aff=false")
                aff = false
            end # If PyPlot is known so we can use it
        end
    end
# ----------------------------------
    println("")
    println("-----------------------------------------")
    println("| time |   primal   |    dual    |  It  |")
    println("-----------------------------------------")
    loop = true
    TiMe = zeros(nbitermax)
    while loop
        tic()
        PDATA.ind += 1

# update of yc from V2
        PDATA.y_v2 = PDATA.yc + PDATA.tau_pwc / rho_y
        proxv2!(PDATA.y_v2, P, W, rho_y, alpha, nb, nw)

# update of yc from phases difference
        PDATA.y_phi = PDATA.yc + PDATA.tau_xic / rho_y
        proxphase!(PDATA.y_phi, Xi, K, rho_y, beta, nb, nw, H, OPTOPT)

# Consensus
        y_tmp = copy(PDATA.yc)
        PDATA.yc = ( PDATA.y_v2 + PDATA.y_phi + PDATA.Fx + (PDATA.tau_xc - (PDATA.tau_pwc + PDATA.tau_xic)) ./ rho_y ) ./ 3

# Object estimation
        x_tmp = copy(PDATA.x)
        PDATA.x,PDATA.Fx = estimx_par(PDATA.x, PDATA.Fx, rho_y, rho_spat, rho_spec, rho_ps,eta ,PDATA.yc, PDATA.z, PDATA.v,
                                      PDATA.w, PDATA.tau_xc, PDATA.tau_s, PDATA.tau_v, PDATA.tau_w, nb, nw, nx, NWvlt,
                                      plan, Wvlt, M, paral)

# update of auxiliary variables
        for n in 1:nw, b in 1:NWvlt
            PDATA.Hx[:, :, n, b] = dwt(PDATA.x[:, :, n], wavelet(Wvlt[b]))
        end

# update of z
        PDATA.z = PDATA.Hx + (PDATA.tau_s / rho_spat)
        PDATA.z = max(1 - ((lambda_spat / rho_spat) ./ abs(PDATA.z)), 0.) .* PDATA.z

# update of v
        tmp = permutedims(rho_spec * PDATA.r - PDATA.tau_r, [3, 1, 2])
        for m in 1:nx, n in 1:nx
            PDATA.Spcdct[m, n, :]= idct(tmp[:, m, n] )
        end
        PDATA.v = ( PDATA.Spcdct + (rho_spec * PDATA.x) + PDATA.tau_v) / (2 * rho_spec)
        vecv = permutedims(PDATA.v, [3, 1, 2])
        for m in 1:nx, n in 1:nx
            PDATA.vHt[m, n, :] = dct(vecv[:, m, n] )
        end

# update of r
        PDATA.r = PDATA.vHt + PDATA.tau_r / rho_spec
        PDATA.r = max(1 - (lambda_spec / rho_spec) ./ abs(PDATA.r), 0) .* PDATA.r

# update of w
        u = PDATA.x + PDATA.tau_w ./ rho_ps
        PDATA.w = max(max(0.0, u) .* mask3D - lambda_L1, 0)

# update of Lagrange multipliers
        PDATA.tau_pwc = PDATA.tau_pwc + rho_y * (PDATA.yc - PDATA.y_v2)
        PDATA.tau_xic = PDATA.tau_xic + rho_y * (PDATA.yc - PDATA.y_phi)
        PDATA.tau_xc = PDATA.tau_xc + rho_y * (PDATA.Fx - PDATA.yc)
        PDATA.tau_s = PDATA.tau_s + rho_spat * (PDATA.Hx - PDATA.z)
        PDATA.tau_v = PDATA.tau_v + rho_spec * (PDATA.x - PDATA.v)
        PDATA.tau_w = PDATA.tau_w + rho_ps * (PDATA.x - PDATA.w)
        PDATA.tau_r = PDATA.tau_r + rho_spec * (PDATA.vHt- PDATA.r)

# stopping criteria
        n1 = norm(vec(PDATA.x -x_tmp))
        n2 = norm(vec(PDATA.yc-y_tmp))
        push!(PDATA.crit1, n1)
        push!(PDATA.crit2, n2)

# Plot and verbose
        if aff&&(PDATA.ind - 1)==(PDATA.count * PDATA.CountPlot)
            OIDATA.PlotFct(PDATA,OIDATA)
            PDATA.count += 1
        end

        if (PDATA.ind >= nbitermax)||( (n1 < eps1)&&(n2 < eps2) )
            loop = false
        end

        @printf("| %02.02f | %02.04e | %02.04e | %04d |\n",toq(), PDATA.crit1[PDATA.ind], PDATA.crit2[PDATA.ind], PDATA.ind)

    end
    return PDATA
end
###################################################################################
# PAINTER MAIN FUNCTION
###################################################################################
function painter(;Folder = "", nbitermax = 1000, nx = 64, lambda_spat = 1/nx^2,
                 lambda_spec = 1/100, lambda_L1 = 0, epsilon = 1e-6,
                 rho_y = 1, rho_spat = 1, rho_spec = 1, rho_ps = 1, alpha = 1,
                 dptype = "all", dpprm = 0,
                 Wvlt  = [WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar],
                 beta = 1, eps1 = 1e-6, eps2 = 1e-6, FOV = 4e-2, mask3D = [], xinit3D = [], indfile = [], indwvl = [],
                 ls = OptimPack.MoreThuenteLineSearch(ftol = 1e-4, gtol = 0.9),
                 scl = OptimPack.SCALING_OREN_SPEDICATO, gat = 1e-6, grt = 1e-6,
                 vt = false, memsize = 100, mxvl = 1000, mxtr = 1000, stpmn = 1e-20,
                 stpmx = 1e+20, PlotFct = painterplotfct, aff = false, CountPlot = 10, admm = true, paral = true)

# Check if mandatory package are installed
    checkPack()
# PAINTER Data Type Creation
    OIDATA = painterinputinit()
    PDATA  = painterdatainit()
    OPTOPT = optiminit(ls, scl, gat, grt, vt, memsize, mxvl, mxtr, stpmn, stpmx)
# PAINTER User parameter validation
    OIDATA = painterinit(OIDATA, Folder, nx, lambda_spat, lambda_spec, lambda_L1, 
                         epsilon, rho_y, rho_spat, rho_spec, rho_ps, alpha, beta,
                         eps1, eps2, FOV, mask3D, xinit3D, Wvlt, paral, 
                         dptype, dpprm, PlotFct)
# OIFITS-FITS Data Read
    println("")
    OIDATA = readoifits(OIDATA, indfile, indwvl)
## PAINTER Matrices creation, Array Initialization
    println("")
    PDATA,OIDATA  = painterarrayinit(PDATA, OIDATA)
# Check, Create PAINTER object and mask initialisation from data or fits
    OIDATA.mask3D   = checkmask(OIDATA.mask3D, OIDATA.nx, OIDATA.nw)
    PDATA.CountPlot = CountPlot
# Main Loop ADMM
    if admm
        PDATA = painteradmm(PDATA, OIDATA, OPTOPT, nbitermax, aff)
    end
    return OIDATA, PDATA, OPTOPT
end
#########################################
function painter(OIDATA::PAINTER_Input,PDATA::PAINTER_Data,OPTOPT::OptOptions,nbitermax::Int,aff::Bool;PlotFct = painterplotfct)
    OIDATA.PlotFct = PlotFct
    PDATA = painteradmm(PDATA, OIDATA, OPTOPT, nbitermax, aff)
    return OIDATA, PDATA, OPTOPT
end
