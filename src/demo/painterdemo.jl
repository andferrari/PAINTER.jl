if(Pkg.installed("PyPlot") == nothing)
    Pkg.add("PyPlot")
end

using PyPlot
close("all")

function plotfunction(PDATA::PAINTER.PAINTER_Data,OIDATA::PAINTER.PAINTER_Input)
    nx = OIDATA.nx
    nw = OIDATA.nw
    wvl = OIDATA.wvl
    FOV = OIDATA.FOV
    x = PDATA.x
    w = PDATA.w .> 0.

    indpix = linspace(-(FOV / 2), (FOV / 2), nx)
    pos = round.(Int, [1, round(Int,nx / 4), round(Int,nx / 2), round(Int,nx * 3 / 4), nx])

    count_y = 0
    count_x = 0
    SubRow  = 6
    SubColumn = 5

    if nw == 30

    for n in 1:nw
        subplot(SubColumn, SubRow, n)
        imshow(x[:, :, n] .* max.(0, w[:, :, n]), origin = "lower")
        titlestring = @sprintf("%2.4f µm", wvl[n] * 1e6)
        title(titlestring)
        xticks([])
        yticks([])
        if( n == (nw + 1 - SubRow + count_x) )
            xticks(collect(pos - 1), round.(Int,indpix[pos] * 100000) / 100)
            xlabel("FOV (mas)")
            count_x += 1
        end
        if(n == (1 + count_y * SubRow))
            yticks(collect(pos - 1), round.(Int,indpix[pos] * 100000) / 100)
            ylabel("FOV (mas)")
            count_y += 1
        end
    end

    else
        NW = ceil(Int, nw/30)
        m = 0
        for n in 1:NW:nw
            m+=1
            subplot(SubColumn, SubRow, m)
            imshow(x[:, :, n] .* max(0, w[:, :, n]), origin = "lower")
            titlestring = @sprintf("%2.4f µm", wvl[n] * 1e6)
            title(titlestring)
            xticks([])
            yticks([])
        if( m == (round(Int,nw/NW)  + 1 - SubRow + count_x) )
            xticks(collect(pos - 1), round(Int,indpix[pos] * 100000) / 100)
            xlabel("FOV (mas)")
            count_x += 1
        end
        if( m == (1 + count_y * SubRow))
            yticks(collect(pos - 1), round(Int,indpix[pos] * 100000) / 100)
            ylabel("FOV (mas)")
            count_y += 1
        end
        end
    end

end

    PlotFct = plotfunction
    FOV = 0.01
    indwvl = 1:6:180

    nx = 64

    eps1 = 1e-4
    eps2 = 1e-4

    admm = true
    aff = true     # plot is enabled

    nbitermax = 1
    savepath = "data.jld"
    dptype = "all"
    dpprm = 5
    Folder = ""
    # first estimate is flat normalized (ones everywhere)
    xinit3D = PAINTER.mask(nx, round(Int, nx/2 ) )
    # initialize algorithm and run admm

    rho_y = 10.
    rho_spat = 4.
    rho_ps = rho_spat
    rho_spec = .5
    lambda_L1 = 5*1e-3
    nbitermax = 1

    OIDATA, PDATA = PAINTER.painter(nbitermax = nbitermax, nx = nx, rho_y = rho_y,
                            rho_spat = rho_spat, lambda_L1 = lambda_L1,
                            rho_spec = rho_spec, rho_ps = rho_ps,
                            eps1 = eps1, eps2 = eps2, FOV = FOV, indwvl = indwvl, admm=admm,
                            PlotFct = PlotFct, aff = aff, dptype = dptype,
                            dpprm = dpprm, Folder = Folder, xinit3D = xinit3D)

# save data struture in .jld files
    println("save data")
    PAINTER.paintersave(savepath,PDATA,OIDATA)

# load data struture in .jld files
    println("load data")
    PDATA, OIDATA = PAINTER.painterload(savepath)

# Warm start of the algorithm to Iteration 100
    nbitermax = 100
    println("PAINTER: Warm Start")
    OIDATA, PDATA = PAINTER.painter(OIDATA,PDATA,nbitermax,true, PlotFct = PlotFct)

    nothing
