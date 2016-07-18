###################################################################################
# Antony Schutz 2015, ANR - POLCA - 2016
# TOOLS for PAINTER
# Schutz, A., Ferrari, A., Mary, D., Soulez, F., Thiébaut, E., Vannier, M.
# PAINTER: a spatiospectral image reconstruction algorithm for optical interferometry
# J. Opt. Soc. Am. A, Vol. 31, N. 11, pp 2334--2345, 2014.
#
# PAINTER Outpt
#
###################################################################################
# using PyPlot
# verb
function painterplotfct(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)
    x = PDATA.x
    if sum(PDATA.w)==0
      w = ones(PDATA.w)
    else
      w = PDATA.w
    end

    crit1 = PDATA.crit1
    crit2 = PDATA.crit2
    eps1 = OIDATA.eps1
    eps2 = OIDATA.eps2
    wvl = OIDATA.wvl
    nx = OIDATA.nx
    nw = OIDATA.nw
    FOV = OIDATA.FOV
    SubColumn,SubRow = createsubplotindex(nw)

    indpix = linspace(-(FOV / 2), (FOV / 2), nx)
    @printf("stop crit, primal: %e (%e) \t dual %e (%e)\n", crit1[end], eps1, crit2[end], eps2)
###################################################################################
# Reconstructed Cube
###################################################################################
    figure(1)
    clf()
    pos     = round(Int, [1, round(Int, nx / 4), round(Int, nx / 2), round(Int, nx * 3 / 4), nx])
    count_y = 0
    count_x = 0
    for n in 1:nw
        subplot(SubColumn, SubRow, n)
        imshow(x[:, :, n] .* max(0., w[:, :, n]), origin ="lower")
        titlestring = @sprintf("%2.4f µm", wvl[n] * 1e6)
        title(titlestring)
        xticks([])
        yticks([])
        if( n == (nw + 1 - SubRow + count_x) )
            xticks(pos - 1, round(Int, indpix[pos] * 100000) / 100)
            xlabel("FOV (mas)")
            count_x += 1
        end
        if(n == (1 + count_y * SubRow))
            yticks(pos - 1, round(Int, indpix[pos] * 100000) / 100)
            ylabel("FOV (mas)")
            count_y += 1
        end
    end
###################################################################################
# primal and dual residuals plot
###################################################################################
    figure(2)
    clf()
    label1 = "Primal Residual"
    label2 = "Dual Residual"
    label3 = "Primal Threshold"
    label4 = "Dual Threshold"
    figure(2)
    plt_1 = plot(crit1, color="blue", label = label1)
    plt_2 = plot(crit2, color="red" , label = label2)
    plt_3 = plot(eps1 * ones(crit1), linestyle = "--", color = "blue", label = label3)
    plt_4 = plot(eps2 * ones(crit2), linestyle = "--", color = "red" , label = label4)
    ax1 = gca()
    ax1[:set_yscale]("log")
    xlabel("Iteration")
    ylabel("Residual (semilog)")
    title( "Primal and Dual residual")
    legend([plt_1, plt_2, plt_3, plt_4],[label1, label2, label3, label4], loc = "upper right", fancybox = "true")
    grid("on")
    #  show()
end
###################################################################################
# PLOT TOOLS
# ---------------------------------------------------------------------------------
# Automatic Subplot index creation
# (SubColumn,SubRow)= SubIndexPlot(nw::Integer)
function createsubplotindex(nw::Integer)
    if(nw == 1)
        SubColumn = 1
        SubRow = 1
    elseif(nw == 2)
        SubColumn = 1
        SubRow = 2
    elseif(nw == 3)
        SubColumn = 1
        SubRow = 3
    elseif(nw == 4)
        SubColumn = 2
        SubRow = 2
    elseif(nw == 5)
        SubColumn = 2
        SubRow = 3
    elseif(nw == 6)
        SubColumn = 2
        SubRow = 3
    elseif(nw == 7)
        SubColumn = 2
        SubRow = 4
    elseif(nw == 8)
        SubColumn = 2
        SubRow = 4
    elseif(nw == 9)
        SubColumn = 3
        SubRow = 3
    elseif(nw == 10)
        SubColumn = 2
        SubRow = 5
    elseif(nw == 11)
        SubColumn = 3
        SubRow = 4
    elseif(nw == 12)
        SubColumn = 3
        SubRow = 4
    elseif(nw == 13)
        SubColumn = 3
        SubRow = 5
    elseif(nw == 14)
        SubColumn = 3
        SubRow = 5
    elseif(nw == 15)
        SubColumn = 3
        SubRow = 5
    elseif(nw == 16)
        SubColumn = 4
        SubRow = 4
    elseif(nw == 18)
        SubColumn = 3
        SubRow = 6
    else
        println(" Automatic subplot, can be")
        SubColumn = round(Int, sqrt(nw))
        SubRow = round(Int, .5+ sqrt(nw))
    end
    return int(SubColumn), int(SubRow)
end
