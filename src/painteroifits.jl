###################################################################################
# Antony Schutz 2015, ANR - POLCA
###################################################################################
# OIFITS RELATED - Read files, create quantities of interest
###################################################################################

# path: read OIFITS data from files included in folder OIDATA.Folder
# OIDATA: structure containing all OIDATA and users informations
# indfile: restrict the number of files (as listed by dir() or ls() )
# ex: indfile = [1 3 5] will only read information from the 1st, the 3rd and the 5th file
# if empty read all the files
# indwvl: restrict the analyzed wavelengths restricted by index, ex:
# indwvl=1:20 --- indwvl= [1 5 7] --- indwvl = 1:5:20
#
# - - - - - - - - - - - - - - - - - - - restreindre les longueurs d'ondes par des valeurs en metre
#
function readoifits(OIDATA::PAINTER_Input,indfile=[],indwvl=[])
# Read Directory
Folder    = OIDATA.Folder
CDir      = readdir(Folder)
# check files fits and oifits
fits      = filter(r"\.fits$", CDir)
oifits    = filter(r"\.oifits$", CDir)
CDir      = vcat(fits,oifits)
NbFiles   = length(CDir)

if NbFiles ==0
  error("No fits or oifits files in folder: $Folder")
end

# files reduction, only files indexed by "indfile" are read
if isempty(indfile)
  indfile = 1:NbFiles
end
OIDATA.indwvl = indwvl
# temporary variable for closure index
LT3           = 0
LFil          = 0

# loop on files
for n         = indfile
# Extract Data from OIFITS or FITS Files
name = string(Folder,"/",CDir[n])
println(name)
master = OIFITS.oifits_load(name)


# # # OI wavelength
nw = 0
for db in OIFITS.oifits_select(master, "OI_WAVELENGTH")
nw = length(oifits_get_eff_wave(db))
end
# Variables initialization
if n == minimum(indfile)
  if isempty(OIDATA.indwvl)==1
    OIDATA.indwvl = 1:nw
  end
  # Initialisation of "empty" array
  OIDATA.Closure_index = Array(Int64,0,3)
  OIDATA.U             = Array(Float64,0,nw)
  OIDATA.V             = Array(Float64,0,nw)
  OIDATA.P             = Array(Float64,0,nw)
  OIDATA.W             = Array(Float64,0,nw)
  OIDATA.T3            = Array(Float64,0,nw)
  OIDATA.T3err         = Array(Float64,0,nw)
  OIDATA.DP            = Array(Float64,0,nw)
  OIDATA.DPerr         = Array(Float64,0,nw)
  OIDATA.wvl           = Array(Float64,0,nw)
end

# # # OI wavelength
for db in OIFITS.oifits_select(master, "OI_WAVELENGTH")
OIDATA.wvl = vcat(OIDATA.wvl,oifits_get_eff_wave(db)')
end

# # # OI VIS
if isempty(OIFITS.oifits_select(master, "OI_VIS"))
error("OI_VIS field is missing, PAINTER needs differential visibilities")
else
for dbvis1 in OIFITS.oifits_select(master, "OI_VIS")
# println("OI_VIS")
OIDATA.DP    = vcat(OIDATA.DP   ,oifits_get_visphi(dbvis1)')
OIDATA.DPerr = vcat(OIDATA.DPerr,oifits_get_visphierr(dbvis1)')
end
end

# # OI VIS 2
if isempty(OIFITS.oifits_select(master, "OI_VIS2"))
error("OI_VIS2 field is missing, PAINTER needs squared visibilities")
else
U_v21         = Array(Float64,0)
for dbvis2 in OIFITS.oifits_select(master,  "OI_VIS2")
U,V   = spatialfrequencies(oifits_get_ucoord(dbvis2),oifits_get_vcoord(dbvis2),oifits_get_eff_wave(dbvis2))
U_v21 = vcat(U_v21,oifits_get_ucoord(dbvis2)) # for closure index
OIDATA.U = vcat(OIDATA.U,U)
OIDATA.V = vcat(OIDATA.V,V)
OIDATA.P = vcat(OIDATA.P,oifits_get_vis2data(dbvis2)')
OIDATA.W = vcat(OIDATA.W,oifits_get_vis2err(dbvis2)')
end
end

# # OI T3
if isempty(OIFITS.oifits_select(master, "OI_T3"))
error("OI_T3 field is missing, PAINTER needs closure phases")
else
U_t31         = Array(Float64,0)
U_t32         = Array(Float64,0)
for dbvis3 in OIFITS.oifits_select(master, "OI_T3")
OIDATA.T3    = vcat(OIDATA.T3   ,oifits_get_t3phi(dbvis3)')
OIDATA.T3err = vcat(OIDATA.T3err,oifits_get_t3phierr(dbvis3)')
# Closure index
U_t31 = oifits_get_u1coord(dbvis3)
U_t32 = oifits_get_u2coord(dbvis3)
(OIDATA.Closure_index, LT3, LFil) = findclosureindex(U_v21,U_t31,U_t32,OIDATA.Closure_index,LT3,LFil)
end
end

end

## FINAL OUTPUT
OIDATA.FilesName = CDir
OIDATA.indfile   = indfile


unitconvert = pi/180
OIDATA.U    = OIDATA.U[:,OIDATA.indwvl]
OIDATA.V    = OIDATA.V[:,OIDATA.indwvl]
OIDATA.P    = OIDATA.P[:,OIDATA.indwvl]
OIDATA.W    = OIDATA.W[:,OIDATA.indwvl]
OIDATA.T3   = OIDATA.T3[:,OIDATA.indwvl]*unitconvert
OIDATA.DP   = OIDATA.DP[:,OIDATA.indwvl]*unitconvert
OIDATA.T3err= OIDATA.T3err[:,OIDATA.indwvl]*unitconvert
OIDATA.DPerr= OIDATA.DPerr[:,OIDATA.indwvl]*unitconvert
OIDATA.wvl  = mean(OIDATA.wvl[:,OIDATA.indwvl],1)
OIDATA.wvl  = vec(OIDATA.wvl)
# # compute Differential Phases as needed for PAINTER and adjust angle unit
OIDATA.nb   = size(OIDATA.U,1)
OIDATA.nw   = size(OIDATA.U,2)
OIDATA.DP   = diffphi(   OIDATA.DP   ,OIDATA.nb,OIDATA.nw)
OIDATA.DPerr= diffphierr(OIDATA.DPerr,OIDATA.nb,OIDATA.nw)
# # Verification des unites - TO BE DONE -- TO BE DONE -- TO BE DONE -- TO BE DONE -- TO BE DONE -- TO BE DONE
# check for angle unit
# DPMat    *= pi/180
# T3Mat    *= pi/180
# K3Mat    *= pi/180
# KPMat    *= pi/180



OIDATA.Xi = vcat(vec(OIDATA.T3)   ,OIDATA.DP)
OIDATA.K  = vcat(vec(OIDATA.T3err),OIDATA.DPerr)
  return OIDATA
end

###################################################################################
# FUNCTIONS FOR OIFITS DATA TRANSFORMATION
###################################################################################
# Spatial Frequencies construction
# ---------------------------------------------------------------------------------
function spatialfrequencies(U_Coord::Vector,V_Coord::Vector,wvl::Vector)
# Create spatiale frequencies from bases coordinates and wavelength
# (U,V) = SpatialFrequencies(U_Coord::Vector,V_Coord::Vector,wvl::Vector)
# U_Coord:  U coordinate of V2 (bases not spatial frequencies) vector (U not divided by wavelength)
# V_Coord:  V coordinate of V2 (bases not spatial frequencies) vector (U not divided by wavelength)
# wvl    :  vector of effective used wavelength
U_3D     = repmat(U_Coord,1,length(wvl))./repmat(wvl',length(U_Coord),1)
V_3D     = repmat(V_Coord,1,length(wvl))./repmat(wvl',length(V_Coord),1)
return U_3D,V_3D
end
# ---------------------------------------------------------------------------------
# PHASE CLOSURE INDEX (sample's index) creation
# ---------------------------------------------------------------------------------
function findclosureindex(V2_U_Coord::Vector,T3_U1_Coord::Vector,T3_U2_Coord::Vector,Closure_index::Matrix,LT3::Integer,LFil::Integer)
# update Closure_index array for each files
# deduce the index in sample related to phases vector, from U1 and U2 coord of phi T3 to U coord of V2
# Closure_index = CreateClosureIndex(V2_U_Coord::Vector,T3_U1_Coord::Vector,T3_U2_Coord::Vector)
# V2_U_Coord:  U coordinate of V2 (bases not spatial frequencies) vector (U not divided by wavelength)
# T3_U1_Coord: U coordinate of T3 first base (not spatial frequencies) vector (U1 not divided by wavelength)
# T3_U2_Coord: U coordinate of T3 second base (not spatial frequencies) vector (U2 not divided by wavelength)
# U3 is given by default
# LT3 and LFil are used to increment correctly information from several files
# Closure_index: index in tabs of associated visibilities
#
utmp   = round(1000*V2_U_Coord )
U1     = round(1000*T3_U1_Coord)
U2     = round(1000*T3_U2_Coord)
U3     = round(1000*(T3_U1_Coord+T3_U2_Coord))
for k  = 1 : length(U1)
  pos1 = LT3+find(x->(x==U1[k]),utmp)
  pos2 = LT3+find(x->(x==U2[k]),utmp)
  pos3 = LT3+find(x->(x==U3[k]),utmp)
  Closure_index = vcat(Closure_index,hcat(pos1,pos2,pos3))
end
LT3    = maximum(Closure_index[:])
LFil   = LFil+length(U1)
return Closure_index, LT3, LFil
end
# ---------------------------------------------------------------------------------
# PAINTER used Differential Visibilities, creation from Differential Visibilities
# ---------------------------------------------------------------------------------
function diffphi(DiffPhi::Array,nb::Int,nw::Int)
# DiffPhi: Differential visibilities as given in OIFITS Data
# DiffPhi is matrix (Nbase*Nwvl)
# All channels reference are assumed to be equal
# defined as in equation 21 of PAINTER (1)
# D1              = DiffPhi[:,1];
# Dr              = DiffPhi[:,2:end];
# D               = repmat(D1,1,size(DiffPhi,2)-1);
# DiffPhiAB       = D - Dr;

HDP             = phasetophasediff(rand(1,1),nw,nb,0,1)
DiffPhiAB       = HDP*vec(DiffPhi)

# println(size(DiffPhiAB))
# println(size(HDP))
# println(size(vec(DiffPhi)))
# println(size(HDP*vec(DiffPhi)))
# println(sum( vec(DiffPhiAB) - HDP*vec(DiffPhi)  ))
end
# ---------------------------
# All channels are assumed to be independent
# variance of sum of random variables = sum of random variables variance
function diffphierr(DiffPhierr::Array,nb::Int,nw::Int)
# D1              = DiffPhierr[:,1]
# Dr              = DiffPhierr[:,2:end]
# D               = repmat(D1,1,size(DiffPhierr,2)-1)
# DiffPhiABerr    = D + Dr

HDP             = phasetophasediff(rand(1,1),nw,nb,0,1)
DiffPhiABerr    = abs(HDP)*vec(DiffPhierr)
# println(sum( vec(DiffPhiABerr) -DiffPhiABerr  ))

end
