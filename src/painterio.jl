###################################################################################
# Antony Schutz 2015, ANR - POLCA
###################################################################################
## Painter Save and Load utilities
function paintersave(savepath::ASCIIString,PDATA::PAINTER_Data,OIDATA::PAINTER_Input,OPTOPT::OptOptions)
    JLD.save(string("OIDATA_",savepath),
        "Folder", OIDATA.Folder, "FilesName", OIDATA.FilesName, "indfile", OIDATA.indfile, "indwvl", OIDATA.indwvl, "wvl", OIDATA.wvl, "U",
        OIDATA.U, "V",OIDATA.V, "P", OIDATA.P, "W", OIDATA.W, "Xi", OIDATA.Xi, "K", OIDATA.K, "Closure_index", OIDATA.Closure_index,
        "nb", OIDATA.nb, "nw",OIDATA.nw, "nx",OIDATA.nx, "FOV", OIDATA.FOV, "lambda_spat", OIDATA.lambda_spat, "lambda_spec",
        OIDATA.lambda_spec, "lambda_L1", OIDATA.lambda_L1, "rho_y", OIDATA.rho_y, "rho_spat", OIDATA.rho_spat, "rho_spec",
        OIDATA.rho_spec, "rho_ps", OIDATA.rho_ps, "alpha", OIDATA.alpha, "beta", OIDATA.beta, "eps1", OIDATA.eps1, "eps2",
        OIDATA.eps2, "epsilon", OIDATA.epsilon, "mask3D", OIDATA.mask3D, "xinit3D", OIDATA.xinit3D, "Wvlt", OIDATA.Wvlt, "paral",
        OIDATA.paral, "T3", OIDATA.T3, "T3err", OIDATA.T3err, "DP", OIDATA.DP, "DPerr", OIDATA.DPerr)#, "PlotFct", OIDATA.PlotFct)

    JLD.save(string("PDATA_",savepath),
        "eta", PDATA.eta,"plan", PDATA.plan, "F3D", PDATA.F3D, "H", PDATA.H, "M", PDATA.M, "x", PDATA.x, "vHt",
        PDATA.vHt,"z", PDATA.z, "Hx", PDATA.Hx, "tau_s", PDATA.tau_s, "w", PDATA.w, "v", PDATA.v, "r", PDATA.r, "tau_w",
        PDATA.tau_w, "tau_v", PDATA.tau_v, "tau_r", PDATA.tau_r, "Spcdct", PDATA.Spcdct, "Fx", PDATA.Fx, "tau_xc",
        PDATA.tau_xc, "tau_pwc", PDATA.tau_pwc, "tau_xic", PDATA.tau_xic,
        "yc", PDATA.yc, "y_v2", PDATA.y_v2, "y_phi", PDATA.y_phi, "crit1", PDATA.crit1, "crit2",
        PDATA.crit2, "ind", PDATA.ind, "CountPlot", PDATA.CountPlot, "count", PDATA.count)#, "nbitermax", PDATA.nbitermax)
end
function painterload(loadpath::ASCIIString)
# OIDATA
    OIDATA = painterinputinit()
    tmp = JLD.load(string("OIDATA_",loadpath))
#    OIDATA.PlotFct = tmp["PlotFct"]
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
    OIDATA.paral = tmp["paral"]

# PDATA
    PDATA  = painterdatainit()
    tmp = JLD.load(string("PDATA_",loadpath))
    PDATA.eta = tmp["eta"]
    PDATA.plan = tmp["plan"]
    PDATA.F3D = tmp["F3D"]
    PDATA.H = tmp["H"]
    PDATA.M = tmp["M"]
    PDATA.x = tmp["x"]
    PDATA.vHt = tmp["vHt"]
    PDATA.z = tmp["z"]
    PDATA.Hx = tmp["Hx"]
    PDATA.tau_s = tmp["tau_s"]
    PDATA.w = tmp["w"]
    PDATA.v = tmp["v"]
    PDATA.r = tmp["r"]
    PDATA.tau_w = tmp["tau_w"]
    PDATA.tau_v = tmp["tau_v"]
    PDATA.tau_r = tmp["tau_r"]
    PDATA.Spcdct = tmp["Spcdct"]
    PDATA.Fx = tmp["Fx"]
    PDATA.tau_xc = tmp["tau_xc"]
    PDATA.tau_pwc = tmp["tau_pwc"]
    PDATA.tau_xic = tmp["tau_xic"]
    PDATA.yc = tmp["yc"]
    PDATA.y_v2 = tmp["y_v2"]
    PDATA.y_phi = tmp["y_phi"]
    PDATA.crit1 = tmp["crit1"]
    PDATA.crit2 = tmp["crit2"]
    PDATA.ind = tmp["ind"]
    PDATA.CountPlot = tmp["CountPlot"]
    PDATA.count = tmp["count"]
# OPTOPT
#   tmp = JLD.load(string("OPTOPT_",loadpath))
#   ls  = tmp["ls"]
#   scl  = tmp["scl"]
#   gat  = tmp["gat"]
#   grt  = tmp["grt"]
#   vt  = tmp["vt"]
#   memsize  = tmp["memsize"]
#   mxvl  = tmp["mxvl"]
#   mxtr  = tmp["mxtr"]
#   stpmn  = tmp["stpmn"]
#   stpmx  = tmp["stpmx"]
#   OPTOPT = optiminit(ls,scl,gat,grt,vt,memsize,mxvl,mxtr,stpmn,stpmx)
    return PDATA, OIDATA
end
