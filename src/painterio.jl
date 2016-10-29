###################################################################################
# Antony Schutz 2015, ANR - POLCA - 2016
###################################################################################
## PAINTER Save and Load utilities
function paintersave(savepath::String,PDATA::PAINTER_Data,OIDATA::PAINTER_Input) #,OPTOPT::OptOptions)
    JLD.save(string("OIDATA_",savepath),
        "Folder", OIDATA.Folder,
        "FilesName", OIDATA.FilesName,
        "indfile", OIDATA.indfile,
        "indwvl", OIDATA.indwvl,
        "wvl", OIDATA.wvl,
        "U", OIDATA.U,
        "V",OIDATA.V,
        "P", OIDATA.P,
        "W", OIDATA.W,
        "Xi", OIDATA.Xi,
        "K", OIDATA.K,
        "Closure_index", OIDATA.Closure_index,
        "nb", OIDATA.nb,
        "nw",OIDATA.nw,
        "nx",OIDATA.nx,
        "FOV", OIDATA.FOV,
        "lambda_spat", OIDATA.lambda_spat,
        "lambda_spec", OIDATA.lambda_spec,
        "lambda_L1", OIDATA.lambda_L1,
        "rho_y", OIDATA.rho_y,
        "rho_y_gamma", OIDATA.rho_y_gamma,
        "rho_y_xi", OIDATA.rho_y_xi,

        "rho_spat", OIDATA.rho_spat,
        "rho_spec", OIDATA.rho_spec,
        "rho_ps", OIDATA.rho_ps,
        "alpha", OIDATA.alpha,
        "beta", OIDATA.beta,
        "eps1", OIDATA.eps1,
        "eps2", OIDATA.eps2,
        "epsilon", OIDATA.epsilon,
        "mask3D", OIDATA.mask3D,
        "xinit3D", OIDATA.xinit3D,
        "Wvlt", OIDATA.Wvlt,
        "T3", OIDATA.T3,
        "T3err", OIDATA.T3err,
        "DP", OIDATA.DP,
        "DPerr", OIDATA.DPerr,
        "dptype", OIDATA.dptype ,
        "dpprm", OIDATA.dpprm,
        "baseNb", OIDATA.baseNb,
        "orderedCluster", OIDATA.orderedCluster,
        "norm", OIDATA.norm
        )

    JLD.save(string("PDATA_",savepath),
        "eta", PDATA.eta,
        # "plan", PDATA.plan,
        "F3D", PDATA.F3D,
        "H", PDATA.H,
        "M", PDATA.M,
        "x", PDATA.x,
        "z", PDATA.z,
        "v", PDATA.v,
        "r", PDATA.r,
        "w", PDATA.w,
        "tau_s", PDATA.tau_s,
        "tau_w", PDATA.tau_w,
        "tau_v", PDATA.tau_v,
        "tau_r", PDATA.tau_r,
        "tau_xc", PDATA.tau_xc,
        "tau_pwc", PDATA.tau_pwc,
        "tau_xic", PDATA.tau_xic,
        "Fx", PDATA.Fx,
        "yc", PDATA.yc,
        "y_v2", PDATA.y_v2,
        "y_phi", PDATA.y_phi,
        "crit1", PDATA.crit1,
        "crit2", PDATA.crit2,
        "ind", PDATA.ind,
        "CountPlot", PDATA.CountPlot,
        "count", PDATA.count
        )
end

function painterload(loadpath::String)
# OIDATA
    OIDATA = painterinputinit()
    tmp = JLD.load(string("OIDATA_",loadpath))
    OIDATA.T3 = tmp["T3"]
    OIDATA.DP = tmp["DP"]
    OIDATA.T3err = tmp["T3err"]
    OIDATA.DPerr = tmp["DPerr"]
    OIDATA.Folder = tmp["Folder"]
    OIDATA.FilesName = tmp["FilesName"]
    OIDATA.indfile = tmp["indfile"]
    OIDATA.indwvl = tmp["indwvl"]
    OIDATA.wvl = tmp["wvl"]
    OIDATA.U = tmp["U"]
    OIDATA.V = tmp["V"]
    OIDATA.P = tmp["P"]
    OIDATA.W = tmp["W"]
    OIDATA.Xi = tmp["Xi"]
    OIDATA.K = tmp["K"]
    OIDATA.Closure_index = tmp["Closure_index"]
    OIDATA.nb = tmp["nb"]
    OIDATA.nw = tmp["nw"]
    OIDATA.nx = tmp["nx"]
    OIDATA.FOV = tmp["FOV"]
    OIDATA.lambda_spat = tmp["lambda_spat"]
    OIDATA.lambda_spec = tmp["lambda_spec"]
    OIDATA.lambda_L1 = tmp["lambda_L1"]

    OIDATA.rho_y = tmp["rho_y"]

    OIDATA.rho_y_gamma = tmp["rho_y_gamma"]
    OIDATA.rho_y_xi = tmp["rho_y_xi"]

    OIDATA.rho_spat = tmp["rho_spat"]
    OIDATA.rho_spec = tmp["rho_spec"]
    OIDATA.rho_ps = tmp["rho_ps"]
    OIDATA.alpha = tmp["alpha"]
    OIDATA.beta = tmp["beta"]
    OIDATA.eps1 = tmp["eps1"]
    OIDATA.eps2 = tmp["eps2"]
    OIDATA.epsilon = tmp["epsilon"]
    OIDATA.mask3D = tmp["mask3D"]
    OIDATA.xinit3D = tmp["xinit3D"]
    OIDATA.Wvlt = tmp["Wvlt"]
    OIDATA.dptype = tmp["dptype"]
    OIDATA.dpprm = tmp["dpprm"]
    OIDATA.baseNb = tmp["baseNb"]
    OIDATA.orderedCluster = tmp["orderedCluster"]
    OIDATA.norm = tmp["norm"]
# PDATA
    DegRd = 2.0 * pi / 360.0
    RadAs = 3600.0 / DegRd
    coef = OIDATA.FOV / RadAs

    PDATA  = painterdatainit()
    tmp = JLD.load(string("PDATA_",loadpath))
    PDATA.eta = tmp["eta"]
    # PDATA.plan = tmp["plan"]
    PDATA.plan = planarray_par(OIDATA.U * coef, OIDATA.V * coef, OIDATA.nx, OIDATA.nw)
    PDATA.F3D = tmp["F3D"]
    PDATA.H = tmp["H"]
    PDATA.M = tmp["M"]
    PDATA.x = tmp["x"]
    PDATA.w = tmp["w"]
    PDATA.tau_s = tmp["tau_s"]
    PDATA.tau_w = tmp["tau_w"]
    PDATA.tau_v = tmp["tau_v"]
    PDATA.tau_r = tmp["tau_r"]
    PDATA.tau_xc = tmp["tau_xc"]
    PDATA.tau_pwc = tmp["tau_pwc"]
    PDATA.tau_xic = tmp["tau_xic"]
    PDATA.Fx = tmp["Fx"]
    PDATA.yc = tmp["yc"]
    PDATA.y_v2 = tmp["y_v2"]
    PDATA.y_phi = tmp["y_phi"]
    PDATA.crit1 = tmp["crit1"]
    PDATA.crit2 = tmp["crit2"]
    PDATA.ind = tmp["ind"]
    PDATA.CountPlot = tmp["CountPlot"]
    PDATA.count = tmp["count"]
    PDATA.z = tmp["z"]
    PDATA.v = tmp["v"]
    PDATA.r = tmp["r"]
    return PDATA, OIDATA
end

# ################################################################
# Added by M. Vannier, 10 Apr. 2015, modified by J. Kluska, may 2015
# ################################################################

function painterfitsexport(savepath::String, PDATA::PAINTER_Data, OIDATA::PAINTER_Input; forceWvlExt = false)
    #  savepath : path + name of ouput file, including ".fits" extension
    # OIDATA :  structure which contains all OIFITS information and user defined parameters.
    # PDATA :  structure which contains output results and diagnostics
    # optional : forceWvlExt (default = false) : boolean to force the writing of a specific "WAVELENGTH" extension.
    #            If "false", that extension is produced only if wavelengths have irregular spacing (within 1% precision).
    # EX:
    #  painterfitsexport("~/my_PAINTER/my_file.fits", PDATA, OIDATA)
    #  painterfitsexport("~/my_PAINTER/my_file.fits", PDATA, OIDATA; forceWvlExt = true)

    # test: write specific wavelength extension or not ?
    # Yes if irregular spacing between wvl or if forceWvlExt = true
    # otherwise: only specify the first wvl and increment in main header
    tol = 0.01    # tolerance = +- 1% on diff(wavelength):
    w = OIDATA.wvl
    write_wvl_ext =  !all( ( diff(w) .< mean(diff(w) * (1.0 + tol)) ) & ( diff(w) .> mean(diff(w)* (1.0 - tol) )) ) | (forceWvlExt)

    f = FITS(savepath, "w")
    write(f, sdata(PDATA.x)) # sdata : SharedArray -> Array

    # to be added in main header:
    fits = f.fitsfile
    in_header = [ [ "CRPIX1", 1, "X ref pixel (R.A.)"],
                [ "CRPIX2", 1, "Y ref pixel (dec)"],
                [ "CRPIX3", 1, "wavelength ref pixel (R.A.)"],
                [ "CDELT1", OIDATA.FOV/OIDATA.nx, "X increment (rad)"],
                [ "CDELT2", OIDATA.FOV/OIDATA.nx, "Y increment (rad)"],
                [ "CDELT3", mean(diff(w)), "wavelength increment (m)"],
                [ "CRVAL1", 0, "X minimum (rad): unknown"],
                [ "CRVAL2", 0, "Y minimum (rad): unknown"],
                [ "CRVAL3", minimum(w), "wavelength minimum (m)"]];
    for i in 1:3:size(in_header)[1]
        # print(bytestring(in_header[i]), " ", float(in_header[i + 1]), " ", bytestring(in_header[i + 2]), "\n")
        fits_write_key(fits, bytestring(in_header[i]), float(in_header[i + 1]), bytestring(in_header[i + 2]))
        # ! Warning comes from  float argument, which should be accepted by fits_write_key method ?!
    end

    # make specific wavelength extension :
    if (write_wvl_ext)
      print(" ...Writing wavelength vector in main header")
      coldefs = [("WAVELENGTH", "1D", "m")];
      fits_create_binary_tbl(fits, size(w)[1], coldefs, "WAVELENGTHS")
close(f)
      fits_write_col(fits::FITSFile, Float64, 1::Integer, 1::Integer, 1::Integer, w)
    end

    # make "INFO" extension (actual content to be discussed)
    coldefs = [ ("mu_spat", "1D", ""), # PDATA.lambda_spat
                ("mu_spec", "1D", ""), # PDATA.lambda_spec
                ("epsilon", "1D", ""), # PDATA.epsilon
                ("crit1", "1D", ""),   # PDATA.crit1
                ("crit2", "1D", "") ]; # PDATA.crit1
    fits_create_binary_tbl(fits, 1, coldefs, "INFO")

    # vector of scalar values for successive columns of "INFO". Should match the description in coldefs
    info_values = [OIDATA.lambda_spat, OIDATA.lambda_spec, OIDATA.epsilon, PDATA.crit1, PDATA.crit2]

    #JKL PART
    for i in 1:size(PDATA.crit1)[1]
    	fits_write_col(fits::FITSFile, 1::Integer, i::Integer, 1::Integer, [info_values[1]])
    	fits_write_col(fits::FITSFile, 2::Integer, i::Integer, 1::Integer, [info_values[2]])
    	fits_write_col(fits::FITSFile, 3::Integer, i::Integer, 1::Integer, [info_values[3]])

    	#fits_write_col(fits::FITSFile, 4::Integer, 1::Integer, 1::Integer, [PDATA.crit1])
    	#fits_write_col(fits::FITSFile, 5::Integer, 1::Integer, 1::Integer, [PDATA.crit2])
    end
    fits_write_col(fits::FITSFile, 4::Integer, 1::Integer, 1::Integer, [PDATA.crit1])
    fits_write_col(fits::FITSFile, 5::Integer, 1::Integer, 1::Integer, [PDATA.crit2])
    #JKL END


    #  for i in 1:size(info_values)[1]
    #	print(i)
    #    fits_write_col(fits::FITSFile, i::Integer, 1::Integer, 1::Integer, [info_values[i]])
    #  end

    # ADDING THE VISIBILITIES IN THE FITS FILE (JKL)
    coldefs2 = [ ("U", "2D", ""), # OIDATA.U
                ("V", "2D", ""), # OIDATA.V
                ("VisRe", "2D", ""), # PDATA.Fx.re
                ("VisIm", "2D", ""), # PDATA.Fx.im
		("Wave","2D", "")];   # OIDATA.wvl
    fits_create_binary_tbl(fits, 1, coldefs2, "VIS")

    for i = 1:OIDATA.nb
	    fits_write_col(fits::FITSFile, 1::Integer, i::Integer, 1::Integer, [OIDATA.U[i,:]])
	    fits_write_col(fits::FITSFile, 2::Integer, i::Integer, 1::Integer, [OIDATA.V[i,:]])
	    fits_write_col(fits::FITSFile, 3::Integer, i::Integer, 1::Integer, [real(PDATA.Fx[i,:])])
	    fits_write_col(fits::FITSFile, 4::Integer, i::Integer, 1::Integer, [imag(PDATA.Fx[i,:])])
	    fits_write_col(fits::FITSFile, 5::Integer, i::Integer, 1::Integer, [OIDATA.wvl])
    end

close(f)
end
