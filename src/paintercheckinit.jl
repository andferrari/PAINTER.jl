###################################################################################
# Antony Schutz 2015, ANR - POLCA - 2016
###################################################################################
function checkPack()

    if( Pkg.installed("OptimPack") == nothing )
        error("OptimPack.jl not installed")
    end

    # NFFT
    if( Pkg.installed("NFFT") == nothing )
        error("NFFT.jl not installed")
    end

    # Wavelets
    if( Pkg.installed("Wavelets") == nothing )
        error("Wavelets.jl not installed")
    end

    # OIFITS
    if( Pkg.installed("OIFITS") == nothing )
        error("HDF5.jl not installed")
    end

    # HDF5
    if( Pkg.installed("HDF5") == nothing )
        error("HDF5.jl not installed")
    end

end
# ---------------------------------------------------------------------------------
# Check 3D mask
# ---------------------------------------------------------------------------------
function checkmask(mask3D::Array,nx::Int,nw::Int)
# check if the mask dimensions are valid
# mask3D: Array of dimensions:
# nx*nx*nw: 3D cube mask given by user
# nx^2 *nw: 3D Matrix
# nx*nx: gray image Matrix
# nx^2: gray image Vector
# empty array leads to default mask: no support contraint only positivity cube of ones
#
# if mask3D is a string it must be the path to a fits file the conversion from fits to data is done in painterinit(...)
#
# case: no mask given
    if isempty(mask3D)
        mask3D = ones(Float64, nx, nx, nw)
    end

    dims = ndims(mask3D)

# case, dims=1 mask3D is gray and must be of size (nx*nx)*1
    if dims == 1
        a = length(mask3D)

        if  a == nx * nx
          mask3D = reshape(repmat(reshape(mask3D[:], nx, nx), 1, nw), nx, nx, nw)

        else
          error(" dims of mask3D =1, mask3D must be nx*nx (nx: Nb pixels)")
        end

# case, dims=3 mask3D is multi color and must be nx*nx*nw
    elseif dims == 3
        (a,b,c) = size(mask3D)

        if(a == nx)&&(b == nx)&&(c == nw)
        else
            error(" dims of mask3D =3, mask3D must be nx*nx*nw (nx: Nb pixels, nw: Nb wavelengths)")
        end

# case, dims=2 mask3D is multi color and must be (nx*nx)*nw
# or it is an image (nx*nx) and must be (nx*nx)*nw
    elseif dims == 2
        (a,b) = size(mask3D)

        if (a == nx)&&(b == nx)
            mask3D = reshape(repmat(mask3D, 1, nw), nx, nx, nw)

        elseif (a == nx*nx)&&(b == nw)
            mask3D = reshape(mask3D, nx, nx, nw)
            println(" Mask3D nx^2 * nw")

        else
            error(" dims of mask3D =2, mask3D is multi color and must be of size (nx*nx)*nw (nx: Nb pixels, nw: Nb wavelengths) or mask3D is an image and must be of size nx*nx (nx: Nb pixels)")
        end

    end
    return mask3D
end
# ---------------------------------------------------------------------------------
# First Estimate
# ---------------------------------------------------------------------------------
# compute the first estimate from xinit, a 3D an array
function inityc(xinit::Array,nb::Int,nx::Int,nw::Int,plan::Array)
# xinit Array given by checkinit(...) several kind of dimensions are allowed or path to a fits file
# nb: number of bases, spatial frequencies per wavelength
# nx: side size of images in pixels
# nw: number of wavelength
# plan: non uniform fft plan
    yc = zeros(Complex{Float64},nb,nw)

    for n in 1:nw
        yc[:,n] = nfft(plan[n], im * 0. + xinit[:, :, n]) / nx
    end

    return yc
end
function checkinit(xinit::Array,nb::Int,nx::Int,nw::Int,plan::Array)
# Check and Compute the first estimate
# xinit: Array of dimensions:
# nx*nx*nw: 3D cube initialization given by user
# nx^2 *nw: 3D Matrix
# nx*nx: gray image Matrix
# nx^2: gray image Vector
# or is complexe visibilities of size nb*nw Matrix or (nb*nx)*1 Vector
#
# nb: number of bases, spatial frequencies per wavelength
# nx: side size of images in pixels
# nw: number of wavelength
# plan: non uniform fft plan
    if isempty(xinit)
        xinit = zeros(nx, nx)
        xinit[round(Int,  (nx / 2) + 1), round(Int, (nx / 2) + 1)] = 1
    end

    if ndims(xinit) == 1
        if size(xinit, 1) == (nx * nx)
            yc = inityc(reshape(repmat(reshape(xinit, nx, nx), 1, nw), nx, nx, nw), nb, nx, nw, plan)

        elseif size(xinit,1) == nb * nw
            yc = reshape(xinit, nb, nw)
        end

    elseif ndims(xinit) == 2

        if(size(xinit, 1) == nx)&&(size(xinit, 2) == nx)
            println("Gray Initial estimate")
            yc = inityc(reshape(repmat(xinit, 1, nw), nx, nx, nw), nb, nx, nw, plan)

        elseif(size(xinit, 1) == (nx * nx))&&(size(xinit, 2) == nw)
            println("Colored vectorized Initial estimate")
            yc = inityc(reshape(xinit, nx, nx, nw), nb, nx, nw, plan)

        elseif(size(xinit, 1) == nb)&&(size(xinit, 2) == nw)
            println("Input is Colored complexe visibilities")
            yc = xinit + im * 0.
        end

    end
    return yc,xinit
end
function checkinit(xinit::ASCIIString,nb::Int,nx::Int,nw::Int,plan::Array)
# method 2 read initial estimate from fits and use method
    xfits = read(FITS(xinit))
    return checkinit(xfits, nb, nx, nw, plan)
end
###################################################################################
function mask(nx::Int,side::Int;choice="square")
# create image 2D of mask for mask3D
# choice = square or disk (disk) default: rect
# nx: size of the constructed image
# side: parameter for mask as radius of disk or half side length or square
    if(side > (nx / 2))
        error("side must be less than nx/2")
    end


    if(choice == "square")
        side2 = (nx / 2) - side
        mask3D = ones(nx, nx)
        mask3D[1:(1 + side2), :] = 0
        mask3D[:, 1:(1 + side2)] = 0
        mask3D[(nx - side2):nx, :] = 0
        mask3D[:, (nx - side2):nx] = 0

    elseif(choice == "disk")
        xg = reshape(repeat([1:nx], outer=[nx]), nx, nx)
        yg = xg'
        rg = sqrt( (xg - .5 - (nx / 2)).^2 + (yg - .5 - (nx / 2)).^2 )
        mask3D = float(abs(rg) .< (side + .5))
    end

    return mask3D
end
###################################################################################
# PAINTER Config
###################################################################################
# Folder        : Folder of fits files
# nx            : size in pixels of image (image of size nx*nx)
# Wvlt          : list of wavelets basis
# lambda__spat  : Spatial Regularization parameter (weight) (Eq.29,31 [1])
# lambda__spec	: Spectral Regularization parameter (weight) (Eq.29,31 [1])
# lambda_L1     : L1 (for star field) constraint, emphasize sparsity of object (Novelty)
# epsilon       : Ridge/tikhonov parameter epilson |X|^2_2 (Eq.29,31 [1])
# rho_y         : ADMM parameter for data fidelity (convergence rate) (Eq.35,50-52 [1])
# rho_spat      : ADMM parameter for Spatial Regularization (convergence rate) (Eq.25,31 [1])
# rho_spec      : ADMM parameter for Spectral Regularization (convergence rate) (Eq.42,55 [1])
# rho_ps        : ADMM parameter for positivity (convergence rate) (Eq.47,54 [1])
# alpha         : Weight for Complexe visibility estimation from V2 (Eq.25,31 [1])
# beta	        : Weight for Complexe visibility estimation from phases difference (Eq.25,31 [1])
# eps1          : Primal Residual stop criterium (3.3 Optimality Conditions and Stopping Criterion , [2])
# eps2          : Dual Residual stop criterium (3.3 Optimality Conditions and Stopping Criterion , [2])
# FOV           : Field Of View User parameter, must be in ArcSecond
# nx            : size of desired output image, in pixel, image is then nx*nx pixels
# mask3D	      : Support constraint: can be a path to a fits file, or to data image (2D or 3D in both case)
# xinit3D       : Initial Estimate: can be a path to a fits file, or to data image (2D or 3D in both case)
###################################################################################
function painterinit(OIDATA::PAINTER_Input,Folder,nx,lambda_spat,lambda_spec,lambda_L1,epsilon,rho_y,rho_spat,rho_spec,rho_ps,alpha,beta,eps1,eps2,FOV,mask3D,xinit3D,Wvlt,paral,dptype,dpprm,PlotFct)
# check if user parameters are valid parameters, correct them if type is not good or replace by default if parameter are not valid
    if(FOV < 0)
        println("FOV must be non negative (default: 0.04 arcsecond)")
        println("initialized to default value")
        FOV = .04
    end

    FOV = FOV + 0.

    if !isinteger(nx)
        println("nx must an integer")
        println("Converted to integer")
    end

    nx = round(Int, nx)

    if(nx < 0)
        println("nx must be non negative (default: 64 pixels)")
        println("initialized to default value")
        nx = 64
    end

    if(lambda_spat < 0)

        println("lambda_spat must be non negative (default: 1.0)")
        println("initialized to default value")
        lambda_spat = 1.
    end

    lambda_spat = lambda_spat + 0.

    if(lambda_spec < 0)
        println("lambda_spec must be non negative (default: 1.0)")
        println("initialized to default value")
        lambda_spec = 1.
    end

    lambda_spec = lambda_spec + 0.

    if(lambda_L1 < 0)
        println("lambda_L1 must be non negative (default: 1e-20)")
        println("initialized to default value")
        lambda_L1 = 1e-20
    end

    lambda_L1 = lambda_L1 + 0.

    if(epsilon < 0)
        println("epsilon must be non negative (default: 1e-6, Ridge Regularization)")
        println("initialized to default value")
        epsilon = 1e-6
    end

    epsilon = epsilon + 0.

    if(rho_y < 0)
        println("rho_y must be non negative (default: 10.0)")
        println("initialized to default value")
        rho_y = 10.
    end

    rho_y = rho_y + 0.

    if(rho_spat < 0)
        println("rho_spat must be non negative (default: 1.0)")
        println("initialized to default value")
        rho_spat = 1.
    end

    rho_spat = rho_spat + 0.

    if(rho_spec < 0)
        println("rho_spec must be non negative (default: 1.0)")
        println("initialized to default value")
        rho_spec = 1.
    end

    rho_spec = rho_spec + 0.

    if(rho_ps < 0)
        println("rho_ps must be non negative (default: 0.5)")
        println("initialized to default value")
        rho_ps = .5
    end

    rho_ps = rho_ps + 0.

    if(alpha < 0)
        println("alpha must be non negative (default: 1.0)")
        println("initialized to default value")
        alpha = 1.
    end

    alpha = alpha + 0.

    if(beta < 0)
        println("beta must be non negative (default: 1.0)")
        println("initialized to default value")
        beta = 1.
    end

    beta = beta + 0.

    if(eps1 < 0)
        println("eps1 must be non negative (default: 1e-4)")
        println("initialized to default value")
        eps1 = 1e-4
    end

    eps1 = eps1 + 0.

    if(eps2 < 0)
        println("eps2 must be non negative (default: 1e-4)")
        println("initialized to default value")
        eps2 = 1e-4
    end

    eps2 = eps2 + 0.

# Check mask3d type
    if typeof(mask3D) == ASCIIString
        mask3dprint = "3D Mask Initialization from fits file: $mask3D"
        f = FITS(mask3D)
        mask3D = read(f[1])

    elseif typeof(mask3D) == Array{Float64, 2}
        Sz = size(mask3D)
        mask3dprint = "3D Mask Initialization from 2D data of size $Sz"

    elseif typeof(mask3D) == Array{Float64, 3}
        Sz = size(mask3D)
        mask3dprint = "3D Mask Initialization from 3D data of size $Sz"

    elseif isempty(mask3D)
        mask3dprint = "3D Mask Initialization from default, no constraint"
    end

# Check xinit3d type
    if typeof(xinit3D) == ASCIIString
        xinit3dprint = "3D Init Initialization from fits file: $xinit3D"
        f = FITS(xinit3D)
        xinit3D = read(f[1])

    elseif typeof(xinit3D) == Array{Float64, 2}
        Sx = size(xinit3D)
        xinit3dprint = "3D Init Initialization from 2D data of size $Sx"

    elseif typeof(xinit3D) == Array{Float64, 3}
        Sx = size(xinit3D)
        xinit3dprint = "3D Init Initialization from 3D data of size $Sx"

    elseif isempty(xinit3D)
        xinit3dprint = "3D Init Initialization from default, centered dirac"
    end

    # Check Differential phases matrix parameters
    if typeof(dptype) == ASCIIString
        if (dptype!="all") && (dptype!="ref") && (dptype!="frame")&& (dptype!="sliding")&& (dptype!="diag")
              println("dptype: all (default), ref, diag, frame, sliding ")
              println("initialized to default value")
              dptype = "all"
        end
    else
        error("dptype must be a string: all (default), ref, frame, sliding, diag ")
    end

    if !isinteger(dpprm)
        println("lambda_ref/horizon must an integer")
        println("Converted to integer")
    end

    dpprm = round(Int,dpprm)

    if(dpprm <= 0)
        println("lambda_ref/horizon must be non negative (default: 1)")
        println("initialized to default value")
        dpprm = 1
    end

# Check data path
    if isempty(Folder)
        tmp = pwd()
        Folder = string(tmp, tmp[1], "OIFITS")
        println(Folder)
    end

    if (typeof(Folder) == ASCIIString)||(typeof(Folder) == UTF8String)
        cpath = Folder

    else
        error("data path is not correct")
    end

# if pwd/OIFITS does not exits so, search in user/.julia/vx.x/PAINTER/src/OIFITS
    if !isdir(cpath)
        pkgfol = Pkg.dir()
        fs = pkgfol[1] # file separator
        cptmp = cpath
        cpath = string(pkgfol, fs, "PAINTER", fs, "src", fs,"OIFITS")
        println("OIFITS files folder does not exists in $cptmp, replaced by default $cpath")
        println("")
    end

    pathprint = cpath

# Check Wavelet input
    if (typeof(Wvlt) <: Array) & all( [typeof(wlt) <: WT.OrthoWaveletClass for wlt in Wvlt])
        Wvltprint = "$Wvlt"
    else
        error("Wavelets list must be a vector of orthogonal wavelets")
    end

    if typeof(PlotFct) != Function
        error("PlotFct must a Function, default plot function: painterplotfct() will be used")
        PlotFct = painterplotfct
    end

    println("OIFits path = $pathprint")
    println("lambda_spat = $lambda_spat ")
    println("lambda_spec = $lambda_spec ")
    println("lambda_L1   = $lambda_L1 ")
    println("epsilon     = $epsilon ")
    println("rho_y       = $rho_y ")
    println("rho_spat    = $rho_spat ")
    println("rho_spec    = $rho_spec ")
    println("rho_ps      = $rho_ps ")
    println("alpha       = $alpha ")
    println("beta	    = $beta ")
    println("eps1        = $eps1 ")
    println("eps2        = $eps2 ")
    println("FOV         = $FOV arcsecond")
    println("nx          = $nx pixels")
    println("mask3D      = $mask3dprint")
    println("xinit3d     = $xinit3dprint")
    println("Wavelets    = $Wvltprint")
    println("Plot Func   = $PlotFct")
    println("dptype      = $dptype")
    println("dpprm       = $dpprm")
    
    OIDATA.PlotFct = PlotFct
    OIDATA.Folder = cpath
    OIDATA.lambda_spat = lambda_spat
    OIDATA.lambda_spec = lambda_spec
    OIDATA.lambda_L1 = lambda_L1
    OIDATA.epsilon = epsilon
    OIDATA.rho_y = rho_y
    OIDATA.rho_spat = rho_spat
    OIDATA.rho_spec = rho_spec
    OIDATA.rho_ps = rho_ps
    OIDATA.alpha = alpha
    OIDATA.beta = beta
    OIDATA.eps1 = eps1
    OIDATA.eps2 = eps2
    OIDATA.FOV = FOV
    OIDATA.nx = nx
    OIDATA.mask3D = mask3D
    OIDATA.xinit3D = xinit3D
    OIDATA.Wvlt = Wvlt
    OIDATA.paral = paral
    OIDATA.dptype = dptype
    OIDATA.dpprm = dpprm
    return OIDATA
end
###################################################################################
# Initialise ADMM Array
function painterarrayinit(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)

    # nx = OIDATA.nx
    # nb = OIDATA.nb
    # nw = OIDATA.nw
    nwvlt = length(OIDATA.Wvlt)
    DegRd = 2.0 * pi / 360.0
    RadAs = 3600.0 / DegRd
    coef = OIDATA.FOV / RadAs

# matrices, constants, plan creation
    PDATA.eta = nwvlt * OIDATA.rho_spat + OIDATA.rho_spec + OIDATA.rho_ps + OIDATA.epsilon
    PDATA.plan = planarray_par(OIDATA.U * coef, OIDATA.V * coef, OIDATA.nx, OIDATA.nw)
    PDATA.F3D = nudft3d_par(OIDATA.U * coef, OIDATA.V * coef, OIDATA.nb, OIDATA.nx, OIDATA.nw)
    PDATA.M  = invmat_par(PDATA.F3D, OIDATA.rho_y, PDATA.eta, OIDATA.nw)
    PDATA.H = phasetophasediff(OIDATA.Closure_index, OIDATA.nw, OIDATA.nb, 1, 1,OIDATA.dptype, OIDATA.dpprm)
# Array Initialization
    PDATA.x = SharedArray(Float64, (OIDATA.nx, OIDATA.nx, OIDATA.nw))
    PDATA.vHt = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw)
    PDATA.z = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw, nwvlt)
    PDATA.Hx = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw, nwvlt)
    PDATA.tau_s = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw, nwvlt)
    PDATA.w = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw)
    PDATA.v = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw)
    PDATA.r = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw)
    PDATA.tau_w = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw)
    PDATA.tau_v = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw)
    PDATA.tau_r = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw)
    PDATA.Spcdct = zeros(Float64, OIDATA.nx, OIDATA.nx, OIDATA.nw)
    PDATA.Fx = SharedArray(Complex{Float64}, (OIDATA.nb, OIDATA.nw))
    PDATA.tau_xc = zeros(Complex128, OIDATA.nb, OIDATA.nw)
    PDATA.tau_pwc = zeros(Complex128, OIDATA.nb, OIDATA.nw)
    PDATA.tau_xic = zeros(Complex128, OIDATA.nb, OIDATA.nw)
    PDATA.y_v2 = zeros(Complex128, OIDATA.nb, OIDATA.nw)
    PDATA.y_phi = zeros(Complex128, OIDATA.nb, OIDATA.nw)
    PDATA.yc, OIDATA.xinit3D = checkinit(OIDATA.xinit3D, OIDATA.nb, OIDATA.nx, OIDATA.nw, PDATA.plan)
    return PDATA,OIDATA
end
