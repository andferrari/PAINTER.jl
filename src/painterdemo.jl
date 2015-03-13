using PyPlot

function myplotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)
        nx = OIDATA.nx
        nw = OIDATA.nw
        wvl = OIDATA.wvl
        FOV = OIDATA.FOV
        x = PDATA.x
        w = PDATA.w.>0.

        indpix = linspace(-FOV/2,FOV/2,nx)
        pos = int([1,round(nx/4),round(nx/2),round(nx*3/4),nx])

        count_y = 0
        count_x = 0
        SubRow  = 6
        SubColumn = 5

        for n =1:nw
                subplot(SubColumn,SubRow,n)
                imshow(x[:,:,n].*max(0,w[:,:,n]), origin ="lower")
                titlestring = @sprintf("%2.4f µm",wvl[n]*1e6)
                title(titlestring)
                xticks([])
                yticks([])
                if( n==(nw+1-SubRow+count_x) )
                        xticks([pos-1],round(indpix[pos]*100000)/100)
                        xlabel("FOV (mas)")
                        count_x+=1
                end
                if(n==(1+count_y*SubRow))
                        yticks([pos-1],round(indpix[pos]*100000)/100)
                        ylabel("FOV (mas)")
                        count_y+=1
                end
        end
end

MyPlotFct     = myplotfunction
MyFOV         = 0.01
Myindwvl      = 1:29
Mynx          = 64
Myeps1        = 1e-4
Myeps2        = 1e-4
Myrho_y       = 10
Myalpha       = 1e4
Mybeta        = 1e5
Myrho_spat    = 4
Myrho_ps      = Myrho_spat
Mylambda_spat = 1e-5
Myrho_spec    = 1/2
Mylambda_spec = 1e-5
Myaff         = true     # plot is enabled
Mynbitermax   = 100
Mypar         = true     # parallel computing is disabled


OIDATA, PDATA, OPTOPT = painter(nbitermax=Mynbitermax, nx=Mynx, lambda_spat=Mylambda_spat=Mylambda_spat, lambda_spec=Mylambda_spec, rho_y= Myrho_y, rho_spat= Myrho_spat, rho_spec= Myrho_spec, rho_ps= Myrho_ps, alpha= Myalpha, beta=Mybeta, eps1=Myeps1, eps2=Myeps2, FOV= MyFOV, indwvl=Myindwvl, ls=OptimPack.MoreThuenteLineSearch(ftol=1e-8,gtol=0.95), scl=OptimPack.SCALING_OREN_SPEDICATO,gat=0,grt=1E-3, vt=false,memsize=100,mxvl=1000,mxtr=1000,stpmn=1E-20,stpmx=1E+20,PlotFct=MyPlotFct,aff=Myaff)



# ###################################################################################
# # Antony Schutz 2015, ANR - POLCA
# ###################################################################################
# # Painter test run on test files
# using PyPlot

# function myplotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)
#   nx = OIDATA.nx
#   nw = OIDATA.nw
#   wvl = OIDATA.wvl
#   FOV = OIDATA.FOV
#   x = PDATA.x
#   w = PDATA.w.>0.

#   indpix = linspace(-FOV/2,FOV/2,nx)
#   pos = int([1,round(nx/4),round(nx/2),round(nx*3/4),nx])

#   count_y = 0
#   count_x = 0
#   SubRow  = 6
#   SubColumn = 5

#   for n =1:nw
#     subplot(SubColumn,SubRow,n)
#     imshow(x[:,:,n].*max(0,w[:,:,n]), origin ="lower")
#     titlestring = @sprintf("%2.4f µm",wvl[n]*1e6)
#     title(titlestring)
#     xticks([])
#     yticks([])
#     if( n==(nw+1-SubRow+count_x) )
#       xticks([pos-1],round(indpix[pos]*100000)/100)
#       xlabel("FOV (mas)")
#       count_x+=1
#     end
#     if(n==(1+count_y*SubRow))
#       yticks([pos-1],round(indpix[pos]*100000)/100)
#       ylabel("FOV (mas)")
#       count_y+=1
#     end
#   end
# end
# # using PyPlot

# iFOV         = 0.01
# iindwvl      = 1:29
# inx          = 64
# ieps1        = 1e-4
# ieps2        = 1e-4
# irho_y       = 10
# iepsilon     = 1e-6
# ialpha       = 1e4
# ibeta        = 1e5
# irho_spat    = 4
# irho_ps      = irho_spat
# ilambda_spat = 1e-5
# irho_spec    = 1/2
# ilambda_L1   = 0
# irho_ss    	 = 1
# ilambda_spec = 1e-5

# iPlotFct     = myplotfunction
# imask3D = mask(inx,int(inx/2 -3))

# iaff                       = true
# inbitermax = 1000
# OIDATA,PDATA,OPTOPT        = painter(nx=inx,lambda_spat=ilambda_spat,lambda_spec=ilambda_spec
#                        ,lambda_L1=ilambda_L1,epsilon=iepsilon,rho_y=irho_y,rho_spat=irho_spat
#                         ,rho_spec=irho_spec,rho_ps=irho_ps,alpha=ialpha,beta=ibeta,eps1=ieps1
#                                  ,eps2=ieps2,FOV=iFOV,mask3D=imask3D,indwvl=iindwvl,aff=iaff,paral=true
#                                  ,ls=OptimPack.MoreThuenteLineSearch(ftol=1e-8,gtol=0.95),
#                                  scl=OptimPack.SCALING_OREN_SPEDICATO,gat=0,grt=1E-3,
#                                  vt=false,memsize=100,mxvl=1000,mxtr=1000,stpmn=1E-20,stpmx=1E+20,PlotFct=iPlotFct);


# # painteradmm(PDATA,OIDATA,OPTOPT,nbitermax=inbitermax,aff-iaff,admm=false)
