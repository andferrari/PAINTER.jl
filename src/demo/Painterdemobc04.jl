using PyPlot

function bc04plotfct(PDATA::PAINTER.PAINTER_Data,OIDATA::PAINTER.PAINTER_Input)
    nx=OIDATA.nx
    nw=OIDATA.nw
    wvl=OIDATA.wvl
    FOV=OIDATA.FOV
    x=PDATA.x
    w=PDATA.w
    if sum(w)==0
      w=1.
    end

    indpix=linspace(-(FOV / 2), (FOV / 2), nx)
    pos=round(Int, [1, round(Int,nx / 4), round(Int,nx / 2), round(Int,nx * 3 / 4), nx])
    X3D= squeeze(sqrt( x .* max(0, w) .* max(x,0) ),3)

    imshow( rotr90(X3D)  , origin = "lower")

    xticks(collect(pos - 1), round(Int,indpix[pos] * 1000) )
    xlabel("FOV (mas)")

    yticks(collect(pos - 1), round(Int,indpix[pos] * 1000) )
    ylabel("FOV (mas)")
end


# To change size of the simulation
# nx pixels

println("")
println(" -----------------------------------------")
println("| 1 Wavelength DEMO Beauty Constest 2004 |")
println(" -----------------------------------------")
println("")

Folder= joinpath(dirname(@__FILE__), "OIFITS_bc04")

savepath="bc04.jld"
CountPlot=25
nbitermax= 2000
PlotFct=bc04plotfct
aff=true
FOV = 0.05 * 1e-3 * 242


# pretty and working, maybe try to adjust a bit, but slowly :)
nx=256

rho_y    = .01
rho_spat = .001
rho_ps=rho_spat

lambda_spat=  .5 * 1e-4
lambda_L1 = 1e-3 * 2

rho_spec=0 # 1 wavelength -> No spectral regularization
epsilon=1e-6 # enough for the pb

## For fast result
xinit3D = PAINTER.mask( nx,round(Int, nx/3 ), choice="disk")
mask3D  = PAINTER.mask( nx,round(Int, nx/2 ))

# initialize algorithm and run admm
OIDATA, PDATA=PAINTER.painter(nbitermax=nbitermax,nx=nx,FOV=FOV,
  rho_y=rho_y, rho_spat=rho_spat,rho_spec=rho_spec,rho_ps=rho_ps,
  lambda_spat=lambda_spat, lambda_L1=lambda_L1,
  xinit3D=xinit3D,mask3D=mask3D,
  PlotFct=PlotFct,aff=aff,CountPlot=CountPlot,Folder=Folder)


# save data struture in .jld files
println("save results of BC2004 data")
PAINTER.paintersave(savepath,PDATA,OIDATA)

# load data struture in .jld files
println("load results of BC2004 data")
PDATA, OIDATA=PAINTER.painterload(savepath)


nothing
