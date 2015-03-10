###################################################################################
# Antony Schutz 2015, ANR - POLCA
###################################################################################
# Painter test run on test files

  include("paintertype.jl")
  include("paintertools.jl")
  include("painterio.jl")
  include("painteroifits.jl")
  include("painterconstmat.jl")
  include("paintercheckinit.jl")
  include("painterplot.jl")

# using PyPlot

iFOV         = 0.01
iindwvl      = 1:29
inx          = 64
ieps1        = 1e-4
ieps2        = 1e-4
irho_y       = 10
iepsilon     = 1e-6
ialpha       = 1e4
ibeta        = 1e5
irho_spat    = 4
irho_ps      = irho_spat
ilambda_spat = 1e-5
irho_spec    = 1/2
ilambda_L1   = 0
irho_ss    	 = 1
ilambda_spec = 1e-5

imask3D = mask(inx,int(inx/2 -3))

iaff                       = true
inbitermax = 1000
OIDATA,PDATA,OPTOPT        = painter(nx=inx,lambda_spat=ilambda_spat,lambda_spec=ilambda_spec
                       ,lambda_L1=ilambda_L1,epsilon=iepsilon,rho_y=irho_y,rho_spat=irho_spat
                        ,rho_spec=irho_spec,rho_ps=irho_ps,alpha=ialpha,beta=ibeta,eps1=ieps1
                                 ,eps2=ieps2,FOV=iFOV,mask3D=imask3D,indwvl=iindwvl,aff=iaff,paral=false);


# painteradmm(PDATA,OIDATA,OPTOPT,nbitermax=inbitermax,aff-iaff,admm=false)
