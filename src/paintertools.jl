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
  # function estimx_par{Tw<:WT.OrthoWaveletClass}(x::SharedArray{Float64,3},Fx::SharedArray{Complex{Float64},2},
    function estimx_par{Tw<:WT.OrthoWaveletClass}(
    rho_y::Float64,rho_spat::Float64,rho_spec::Float64,rho_ps::Float64,eta::Float64,
    yc::Array{Complex{Float64},2},z::Array{Float64,4},v::Array{Float64,3},w::Array{Float64,3},
    tau_xc::Array{Complex{Float64},2},tau_s::Array{Float64,4},tau_v::Array{Float64,3},tau_w::Array{Float64,3},
    nb::Int,nw::Int,nx::Int,NWvlt::Int,plan::Array{NFFT.NFFTPlan{2,0,Float64},1},Wvlt::Array{Tw,1},M::Array{Array{Complex{Float64},2},1})
  # Estimate the constrained, regularized 3D images from complexe visibilities
  # step IV of PAINTER [0,1]
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
    matv = 0
    matw = 0
  # parallel sum of wavelet basis
    wvd = SharedArray{Float64}(nx,nx,nw,NWvlt)
    @sync @parallel for ind in 1:nw*NWvlt
        n,b = ind2sub((nw,NWvlt),ind)
        wvd[:,:,n,b] = idwt(matz[:, :, n, b] , wavelet(Wvlt[b]) )
    end
    matz = 0
    wvd = sum(wvd,4)
  # parallel image reconstruction
    x = SharedArray{Float64}(nx,nx,nw)
    Fx = SharedArray{Complex{Float64}}(nb, nw)
    @sync @parallel for n in 1:nw
        xtmp = (nfft_adjoint(plan[n], maty[:, n]) / nx) + Reg[:, :, n] + wvd[:, :, n] #+ wvd[n]
        xfst = nfft_adjoint(plan[n], M[n] * (nfft(plan[n], xtmp) / nx)) / nx
        xtmp = (xtmp - (rho_y / eta) * xfst) / eta
        x[:,:,n] = real(xtmp)
        Fx[:,n]  = nfft(plan[n], xtmp) / nx
    end
    return x, Fx
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
  function proxv2!(y_v2::Array,P::Array,W::Array,rho_y::Real,alpha::Real,nb::Int,nw::Int)
    mod_y = convert(SharedArray,abs.(y_v2))
    ang_y = angle.(y_v2)
    tmp1 = W.*rho_y /(4 * alpha)
    tmp2 = alpha ./ W
    @sync @parallel for z in 1:nb*nw
        m,n = ind2sub((nb,nw),z)
        sol = max.(0., paintercubicroots( tmp1[m,n] - P[m,n], tmp1[m,n] .* mod_y[m,n]))
        cst = tmp2[m,n] .* (P[m,n] - sol.^2).^2 + .5 .* rho_y * (sol - mod_y[m,n] ).^2
        (a,b) = findmin(cst)
        mod_y[m,n] = sol[b]
    end
    y_v2[:] = mod_y .* exp.(im .* ang_y)
  end
  # ---------------------------------------------------------------------------------
  # ----- Cardano's formula
  function paintercubicroots(c::Real,d::Real)
  #-----   special  case  in  Painter -----
  # a = 1.
  # b = 0.
  # d = -d
  #-----  -----  -----  -----  -----  -----
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
    q = c / 3.
    r = d / 2.

    discriminant = q^3 + r^2
    if discriminant >= 0        # We have 1 real root and 2 imaginary
        s = realcuberoot(r + sqrt(discriminant))
        t = realcuberoot(r - sqrt(discriminant))
        root = s + t     # Just calculate the real root
    else                        # We have 3 real roots
        # In this case (r + sqrt(discriminate)) is complex so the following
        # code constructs the cube root of this complex quantity
        rho = sqrt(r^2 - discriminant)
        cubeRootrho = realcuberoot(rho)    # Cube root of complex magnitude
        thetaOn3 = acos(r / rho) / 3       # Complex angle/3
        crRhoCosThetaOn3 = cubeRootrho * cos(thetaOn3)
        crRhoSinThetaOn3 = cubeRootrho * sin(thetaOn3)
        root = zeros(3)
        root[1] = 2 * crRhoCosThetaOn3
        root[2] = -crRhoCosThetaOn3 - sqrt(3) * crRhoSinThetaOn3
        root[3] = -crRhoCosThetaOn3 + sqrt(3) * crRhoSinThetaOn3
    end
  return root
  end
  # ---------------------------------------------------------------------------------
  # Proximal operator for phases difference (Section 5.2 PAINTER)
  # for parallel
  function proxphase(y_phi::Matrix,Xi::Vector,K::Vector,rho_y::Real,beta::Real
                  ,nb::Int,nw::Int,H::SparseMatrixCSC)
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
    gam_t = abs.(y_t)
    phi_t = angle.(y_t)
    phi_0 = angle.(y_t)

    function cost!{T<:Real}(x_phi::Array{T,1}, g_phi::Array{T,1})
        return costgradphi!(x_phi, g_phi, gam_t, phi_t, y_t, Xi, K, beta, rho_y, H)
    end

    phi = OptimPack.vmlmb(cost!, phi_0, mem = memsize, verb = vt
                          , gatol = gat, grtol = grt, maxeval = mxvl
                          , maxiter = mxtr, lnsrch = ls)

    if phi!=nothing
        Ek = phi_t - phi
        gam = max.(0.0, gam_t .* cos.(Ek))
        y_phi[:] = reshape(gam .* exp.(im * phi), nb, nw)
    else
        y_phi[:] = reshape(gam_t .* exp.(im * phi_t), nb, nw)
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
    gam = max.(0.0, gam_t .* cos.(Ek))
    dphi = H * x_phi - Xi;
    yest = gam .* exp.(im * x_phi)
    w1 = -sum(K .* cos.(dphi))
    w2 = sum(abs.(yest - y_t).^2)
    g_phi[:]= beta .* H' * (sin.(dphi) .* K) - rho_y .* gam .* gam_t .* sin.(Ek)
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
  function painteradmm(PDATA::PAINTER_Data,OIDATA::PAINTER_Input,nbitermax::Int,aff::Bool)
    const nx = OIDATA.nx
    const nb = OIDATA.nb
    const nw = OIDATA.nw
    const wvl = OIDATA.wvl
    const lambda_spat = OIDATA.lambda_spat
    const lambda_spec = OIDATA.lambda_spec
    const lambda_L1 = OIDATA.lambda_L1
    const epsilon = OIDATA.epsilon
    const rho_y = OIDATA.rho_y

    const rho_y_gamma = OIDATA.rho_y_gamma
    const rho_y_xi = OIDATA.rho_y_xi

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
    const eta = PDATA.eta
    const plan = PDATA.plan
    const F3D = PDATA.F3D
    const M = PDATA.M
    const NWvlt = length(Wvlt)
    const H = PDATA.H
    const baseNb = OIDATA.baseNb
    const Nt3indep = length(OIDATA.baseNb)
    yphidict = SharedArray{Complex128}(nb,nw)
    Spcdct = convert( SharedArray, zeros( nx, nx, nw))
    vHt = convert( SharedArray, zeros(nx, nx, nw))
    Hx = convert( SharedArray, zeros( nx, nx, nw, NWvlt))

  # ----------------------------------
  # Check if Pyplot is used to graphics
    if aff
  # check if Pyplot is installed
        if( Pkg.installed("PyPlot") == nothing )
            println("")
            println("PyPlot is not installed: Pkg.add(''PyPlot''), aff=false ")
            aff = false
        end
    end
  # ----------------------------------
        println(" ")
        println("------------------------ ")
        println(" Nb phases Cluster: " , Nt3indep)
        println("------------------------ ")
        println("VMLM will run with file: ")
        println("------------------------ ")
        println("memsize: ", memsize )
        println("verb: ", vt )
        println("grtol: ", grt )
        println("gatol: ", gat )
        println("maxeval: ", mxvl )
        println("maxiter: ", mxtr )
        #println("stpmin: ", stpmn )
        #println("stpmax: ", stpmx )
        #println("scaling: ", scl )
        println("lnsrch: ", ls )
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
        PDATA.y_v2 = PDATA.yc + PDATA.tau_pwc / rho_y_gamma
        proxv2!(PDATA.y_v2, P, W, rho_y_gamma, alpha, nb, nw)

  # update of yc from phases difference
        y_phidat = PDATA.yc + PDATA.tau_xic / rho_y_xi
        @sync @parallel for n in 1:Nt3indep
            yphidict[baseNb[n],:] = proxphase(y_phidat[baseNb[n],:], Xi[n], K[n]
                     , rho_y_xi, beta, length( baseNb[n] ), nw, H[n] )
        end
        PDATA.y_phi = yphidict
        y_phidat = 0
  # Consensus
        y_tmp = copy(PDATA.yc)
        PDATA.yc = ( PDATA.y_v2 - PDATA.tau_xic./rho_y_xi + PDATA.y_phi - PDATA.tau_pwc./rho_y_gamma  + PDATA.Fx + PDATA.tau_xc./ rho_y   ) ./ 3

  # Object estimation
        x_tmp = copy(PDATA.x)
        PDATA.x,PDATA.Fx = estimx_par(rho_y, rho_spat, rho_spec, rho_ps,eta ,PDATA.yc, PDATA.z, PDATA.v, PDATA.w,
                                      PDATA.tau_xc, PDATA.tau_s, PDATA.tau_v, PDATA.tau_w, nb, nw, nx, NWvlt,
                                      plan, Wvlt, M)
  # Cube of flux
        flux_cube = 1. # vec(sum(OIDATA.P,1))

  # update of auxiliary variables
        if rho_spat >0
            tmpx = copy(PDATA.x)
            @sync @parallel for ind in 1:nw*NWvlt
                n,b = ind2sub((nw,NWvlt),ind)
                Hx[:, :, n, b] = dwt(tmpx[:, :, n], wavelet(Wvlt[b]))
            end
            tmpx = 0
  # update of z
            PDATA.z = Hx + (PDATA.tau_s / rho_spat)
            PDATA.z = max.(1 - ((lambda_spat / rho_spat) ./ abs.(PDATA.z)), 0.) .* PDATA.z

            # for n in 1:OIDATA.nw
            #     PDATA.z[:,:,n,:] = max(1 - ((flux_cube[n].*lambda_spat / rho_spat) ./ abs(PDATA.z[:,:,n,:])), 0.) .* PDATA.z[:,:,n,:]
            # end
        end
  # update of v
        if rho_spec >0
            tmp = permutedims(rho_spec * PDATA.r - PDATA.tau_r, [3, 1, 2])
            @sync @parallel for ind in 1:nx*nx
                m,n = ind2sub((nx,nx),ind)
                Spcdct[m, n, :]= idct(tmp[:, m, n] )
            end
            tmp = 0

            PDATA.v = ( Spcdct + (rho_spec * PDATA.x) + PDATA.tau_v) / (2 * rho_spec)
            vecv = permutedims(PDATA.v, [3, 1, 2])
            @sync @parallel for ind in 1:nx*nx
                m,n = ind2sub((nx,nx),ind)
                vHt[m, n, :] = dct(vecv[:, m, n] )
            end
            vecv = 0
  # update of r
            PDATA.r = vHt + PDATA.tau_r / rho_spec
            PDATA.r = max.(1 - (lambda_spec / rho_spec) ./ abs.(PDATA.r), 0) .* PDATA.r
        end
  # update of w
        if rho_ps>0
            u = PDATA.x + PDATA.tau_w ./ rho_ps
            PDATA.w = max.(max.(0.0,u) .* mask3D - lambda_L1, 0)
            # for n in 1:OIDATA.nw
            #     PDATA.w[:,:,n] = max(max(0.0,u[:,:,n]) .* mask3D[:,:,n] - flux_cube[n].*lambda_L1, 0)
            # end
            PDATA.tau_w = PDATA.tau_w + rho_ps * (PDATA.x - PDATA.w)
        end
  # update of Lagrange multipliers
        PDATA.tau_pwc = PDATA.tau_pwc + rho_y_gamma * (PDATA.yc - PDATA.y_v2)
        PDATA.tau_xic = PDATA.tau_xic + rho_y_xi * (PDATA.yc - PDATA.y_phi)
        PDATA.tau_xc = PDATA.tau_xc + rho_y * (PDATA.Fx - PDATA.yc)
        if rho_spat >0
            PDATA.tau_s = PDATA.tau_s + rho_spat * (Hx - PDATA.z)
        end
        if rho_spec >0
            PDATA.tau_v = PDATA.tau_v + rho_spec * (PDATA.x - PDATA.v)
            PDATA.tau_r = PDATA.tau_r + rho_spec * (vHt- PDATA.r)
        end

  # stopping criteria
        n1 = norm(vec(PDATA.x -x_tmp))
        n2 = norm(vec(PDATA.yc-y_tmp))
        push!(PDATA.crit1, n1)
        push!(PDATA.crit2, n2)

  # Plot and verbose
        if aff&&( (PDATA.ind )==1 || (PDATA.ind )==( (PDATA.count-1) * PDATA.CountPlot) )
            OIDATA.PlotFct(PDATA,OIDATA)
            PDATA.count += 1
        end
        if (PDATA.ind >= nbitermax)||( (n1 < eps1)&&(n2 < eps2) )
            loop = false
        end
        @printf("| %02.02f | %02.04e | %02.04e | %04d |\n",toq(), PDATA.crit1[PDATA.ind], PDATA.crit2[PDATA.ind], PDATA.ind)
    end
    PDATA.x = copy(PDATA.x)
    PDATA.Fx = copy(PDATA.Fx)
    PDATA.y_v2 = copy(PDATA.y_v2)
    PDATA.y_phi = copy(PDATA.y_phi)
    return PDATA

  end
  ###################################################################################
  # PAINTER MAIN FUNCTION
  ###################################################################################
  function painter(;Folder = "", nbitermax = 1000, nx = 64, lambda_spat = 1e-3,
                 lambda_spec = 1e-3, lambda_L1 = 1e-3, epsilon = 1e-6,
                 rho_y = 10., rho_spat = 1., rho_spec = 1., rho_ps = 1., alpha = 1., beta = 1.,
                 dptype = "all", dpprm = 0, rho_y_gamma = 10., rho_y_xi = 10.,
                 Wvlt  = [WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar],
                 eps1 = 1e-6, eps2 = 1e-6, FOV = 4e-2, mask3D = [], xinit3D = [], indfile = [], indwvl = [],
                 PlotFct = painterplotfct, aff = false, CountPlot = 10, admm = true, flux = 0, autoinit="true")
  # Check if mandatory package are installed
    checkPack()
  # PAINTER Data Type Creation
    OIDATA = painterinputinit()
    PDATA  = painterdatainit()
  # PAINTER User parameter validation
    OIDATA = painterinit(OIDATA, Folder, nx, lambda_spat, lambda_spec, lambda_L1, 
                         epsilon, rho_y, rho_y_gamma, rho_y_xi, rho_spat, rho_spec, rho_ps, alpha, beta,
                         eps1, eps2, FOV, mask3D, xinit3D, Wvlt,  
                         dptype, dpprm, indwvl, PlotFct)

  # OIFITS-FITS Data Read
    println("")
    OIDATA = readoifits(OIDATA, indfile, indwvl)
  ## PAINTER Matrices creation, Array Initialization
    println("")
    PDATA,OIDATA = painterarrayinit(PDATA, OIDATA)
  # Check, Create PAINTER object and mask initialisation from data or fits
    OIDATA.mask3D = checkmask(OIDATA.mask3D, OIDATA.nx, OIDATA.nw)
  # Initialise Data and xinit from V2 and flux
    if autoinit == "true"
        PDATA,OIDATA = painterautoparametersinit(PDATA,OIDATA)
    end
  # initialise Lagrange multipliers for warm start from initial estimates
    PDATA,OIDATA = painterlagrangemultipliersinit(PDATA,OIDATA)

    PDATA.CountPlot = CountPlot

    # Main Loop ADMM
    if admm
        # PDATA = painteradmm(PDATA, OIDATA, OPTOPT, nbitermax, aff)
        PDATA = painteradmm(PDATA, OIDATA, nbitermax, aff)
    end
    return OIDATA, PDATA
  end
  #########################################
  function painter(OIDATA::PAINTER_Input,PDATA::PAINTER_Data,nbitermax::Int,aff::Bool;PlotFct = painterplotfct)
    OIDATA.PlotFct = PlotFct
    PDATA = painteradmm(PDATA, OIDATA, nbitermax, aff)
    return OIDATA, PDATA
  end

  #########################################
  # DEMO
  function painterdemo()
    painterdemo(0)
  end

  function painterdemo(demo)
    if demo == 0
        include( joinpath(dirname(@__FILE__), "demo", "painterdemo.jl") )
    elseif (demo==1)||(demo=="gravity")
        include( joinpath(dirname(@__FILE__), "demo", "painterdemogravity.jl") )
    elseif (demo==3)||(demo=="matiss")
        include( joinpath(dirname(@__FILE__), "demo", "painterdemomatiss.jl") )
      elseif (demo==2)||(demo=="bc04")
          include( joinpath(dirname(@__FILE__), "demo", "painterdemobc04.jl") )
    end
  end
