using PyPlot

function plotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)
    nx = OIDATA.nx
    nw = OIDATA.nw
    wvl = OIDATA.wvl
    FOV = OIDATA.FOV
    x = PDATA.x
    w = PDATA.w .> 0.

    indpix = linspace(-(FOV / 2), (FOV / 2), nx)
    pos = round(Int, [1, round(Int,nx / 4), round(Int,nx / 2), round(Int,nx * 3 / 4), nx])

    count_y = 0
    count_x = 0
    SubRow  = 6
    SubColumn = 5

    for n in 1:nw
        subplot(SubColumn, SubRow, n)
        imshow(x[:, :, n] .* max(0, w[:, :, n]), origin = "lower")
        titlestring = @sprintf("%2.4f µm", wvl[n] * 1e6)
        title(titlestring)
        xticks([])
        yticks([])
        if( n == (nw + 1 - SubRow + count_x) )
            xticks(collect(pos - 1), round(Int,indpix[pos] * 100000) / 100)
            xlabel("FOV (mas)")
            count_x += 1
        end
        if(n == (1 + count_y * SubRow))
            yticks(collect(pos - 1), round(Int,indpix[pos] * 100000) / 100)
            ylabel("FOV (mas)")
            count_y += 1
        end
    end
end

    PlotFct = plotfunction
    FOV = 0.01
    indwvl = 1:30
    nx = 64
    eps1 = 1e-4
    eps2 = 1e-4
    rho_y = 10.
    alpha = 1e4
    beta = 1e5
    rho_spat = 4.
    rho_ps = rho_spat
    lambda_spat = 1e-5
    rho_spec = 0.5
    lambda_spec = 1e-5
    aff = true     # plot is enabled
    nbitermax = 1
    paral = true     # parallel computing is enabled
    savepath = "data.jld"
    dptype = "sliding"
    dpprm = 5
    Folder = ""
# initialize algorithm and run admm
    OIDATA, PDATA, OPTOPT = painter(nbitermax = nbitermax, nx = nx, lambda_spat = lambda_spat,
                                lambda_spec = lambda_spec, rho_y = rho_y, rho_spat = rho_spat,
                                rho_spec = rho_spec, rho_ps = rho_ps, alpha = alpha, beta = beta,
                                eps1 = eps1, eps2 = eps2, FOV = FOV, indwvl = indwvl,
                                ls = OptimPack.MoreThuenteLineSearch(ftol = 1e-8, gtol = 0.95),
                                scl = OptimPack.SCALING_OREN_SPEDICATO, gat = 0, grt = 1e-3,
                                vt = false, memsize = 100, mxvl = 1000, mxtr = 1000, stpmn = 1e-20,
                                stpmx = 1e+20, PlotFct = PlotFct, aff = aff, paral = paral,
                                dptype = dptype, dpprm = dpprm, Folder = Folder)

# save data struture in .jld files
    paintersave(savepath,PDATA,OIDATA,OPTOPT)

# load data struture in .jld files
    PDATA, OIDATA = painterload(savepath)

# Warm start of the algorithm
    OIDATA, PDATA, OPTOPT = painter(OIDATA,PDATA,OPTOPT,100,true, PlotFct = PlotFct)
