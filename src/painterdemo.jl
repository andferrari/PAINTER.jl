using PyPlot

function myplotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)
    nx = OIDATA.nx
    nw = OIDATA.nw
    wvl = OIDATA.wvl
    FOV = OIDATA.FOV
    x = PDATA.x
    w = PDATA.w .> 0.

    indpix = linspace(-(FOV / 2), (FOV / 2), nx)
    pos = int([1, round(nx / 4), round(nx / 2), round(nx * 3 / 4), nx])

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
            xticks([pos - 1], round(indpix[pos] * 100000) / 100)
            xlabel("FOV (mas)")
            count_x += 1
        end
        if(n == (1 + count_y * SubRow))
            yticks([pos - 1], round(indpix[pos] * 100000) / 100)
            ylabel("FOV (mas)")
            count_y += 1
        end
    end
end

    MyPlotFct = myplotfunction
    MyFOV = 0.01
    Myindwvl = 1:30
    Mynx = 64
    Myeps1 = 1e-4
    Myeps2 = 1e-4
    Myrho_y = 10.
    Myalpha = 1e4
    Mybeta = 1e5
    Myrho_spat = 4.
    Myrho_ps = Myrho_spat
    Mylambda_spat = 1e-5
    Myrho_spec = 0.5
    Mylambda_spec = 1e-5
    Myaff = true     # plot is enabled
    Mynbitermax = 1
    Mypar = true     # parallel computing is enabled
    savepath = "mydata.jld"

# initialize algorithm and run admm
    OIDATA, PDATA, OPTOPT = painter(nbitermax = Mynbitermax, nx = Mynx, lambda_spat = Mylambda_spat,
                                lambda_spec = Mylambda_spec, rho_y = Myrho_y, rho_spat = Myrho_spat,
                                rho_spec = Myrho_spec, rho_ps = Myrho_ps, alpha = Myalpha, beta = Mybeta,
                                eps1 = Myeps1, eps2 = Myeps2, FOV = MyFOV, indwvl = Myindwvl,
                                ls = OptimPack.MoreThuenteLineSearch(ftol = 1e-8, gtol = 0.95),
                                scl = OptimPack.SCALING_OREN_SPEDICATO, gat = 0, grt = 1e-3,
                                vt = false, memsize = 100, mxvl = 1000, mxtr = 1000, stpmn = 1e-20,
                                stpmx = 1e+20, PlotFct = MyPlotFct, aff = Myaff)

# save data struture in .jld files
    paintersave(savepath,PDATA,OIDATA,OPTOPT)

# load data struture in .jld files
    PDATA, OIDATA = painterload(savepath)

# Warm start of the algorithm
    OIDATA, PDATA, OPTOPT = painter(OIDATA,PDATA,OPTOPT,100,true, PlotFct = MyPlotFct)
