###################################################################################
# Antony Schutz 2015, ANR - POLCA - 2016
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

# Check for DP and T3 must be in data at least
    OIDATA.isDP = !(length(indwvl)==1)
    indwvl = collect(indwvl)
# Read Directory
    Folder    = OIDATA.Folder
    CDir      = readdir(Folder)

# check files fits and oifits
    fits      = filter(r"\.fits$", CDir)
    oifits    = filter(r"\.oifits$", CDir)
    CDir      = vcat(fits,oifits)
    NbFiles   = length(CDir)

    if NbFiles == 0
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
# Extract Data from OIFITS or FITS Files

    for n in indfile
        name = string(Folder, "/", CDir[n])
        println(name)
        master = OIFITS.load(name)
# # # OI wavelength
        nw = 0
        for db in OIFITS.select(master, "OI_WAVELENGTH")
            nw = length(OIFITS.get_eff_wave(db))
        end
# Variables initialization
        if n == minimum(indfile)
            if isempty(OIDATA.indwvl)
                OIDATA.indwvl = 1:nw
            end
# Initialisation of "empty" array
            OIDATA.Closure_index = Array{Int}(0, 3)
            OIDATA.U = Array{Float64}(0, nw)
            OIDATA.V = Array{Float64}(0, nw)
            OIDATA.P = Array{Float64}(0, nw)
            OIDATA.W = Array{Float64}(0, nw)
            OIDATA.T3 = Array{Float64}(0, nw)
            OIDATA.T3err = Array{Float64}(0, nw)
            OIDATA.DP = Array{Float64}(0, nw)
            OIDATA.DPerr = Array{Float64}(0, nw)
            OIDATA.wvl = Array{Float64}(0, nw)
        end

# # # OI wavelength
        for db in OIFITS.select(master, "OI_WAVELENGTH")
            OIDATA.wvl = vcat(OIDATA.wvl,OIFITS.get_eff_wave(db)')
        end
# # # OI VIS
        if OIDATA.isDP == 1
            if isempty(OIFITS.select(master, "OI_VIS"))
                println("OI_VIS field is missing, PAINTER needs differential visibilities")
                println(" --- will try to work only with closure phases ")
                OIDATA.isDP = 0
            else

                for dbvis1 in OIFITS.select(master, "OI_VIS")
                    OIDATA.DP = vcat(OIDATA.DP, OIFITS.get_visphi(dbvis1)')
                    OIDATA.DPerr = vcat(OIDATA.DPerr, OIFITS.get_visphierr(dbvis1)')
                end
            end
        end

# # OI VIS 2
        if isempty(OIFITS.select(master, "OI_VIS2"))
            error("OI_VIS2 field is missing, PAINTER needs squared visibilities")

        else
            U_v21 = Array{Float64}(0)

            for dbvis2 in OIFITS.select(master,  "OI_VIS2")
                U,V  = spatialfrequencies(OIFITS.get_ucoord(dbvis2), OIFITS.get_vcoord(dbvis2), OIFITS.get_eff_wave(dbvis2))
                U_v21 = vcat(U_v21, OIFITS.get_ucoord(dbvis2)) # for closure index
                OIDATA.U = vcat(OIDATA.U, U)
                OIDATA.V = vcat(OIDATA.V, V)
                OIDATA.P = vcat(OIDATA.P, OIFITS.get_vis2data(dbvis2)')
                OIDATA.W = vcat(OIDATA.W, OIFITS.get_vis2err(dbvis2)')
            end
        end

# # OI T3
        if isempty(OIFITS.select(master, "OI_T3"))
            error("OI_T3 field is missing, PAINTER needs closure phases")

        else
            U_t31 = Array{Float64}(0)
            U_t32 = Array{Float64}(0)
            for dbvis3 in OIFITS.select(master, "OI_T3")
                OIDATA.T3 = vcat(OIDATA.T3, OIFITS.get_t3phi(dbvis3)')
                OIDATA.T3err = vcat(OIDATA.T3err, OIFITS.get_t3phierr(dbvis3)')
# Closure index
                U_t31 = OIFITS.get_u1coord(dbvis3)
                U_t32 = OIFITS.get_u2coord(dbvis3)
                (OIDATA.Closure_index, LT3, LFil) = findclosureindex(U_v21, U_t31, U_t32, OIDATA.Closure_index, LT3, LFil)
            end
        end
    end

## FINAL OUTPUT - Wavelength reduction
    unitconvert = pi / 180
    OIDATA.FilesName = CDir
    OIDATA.indfile = indfile
    OIDATA.U = OIDATA.U[:, OIDATA.indwvl]
    OIDATA.V = OIDATA.V[:, OIDATA.indwvl]
    OIDATA.P = OIDATA.P[:, OIDATA.indwvl]
    OIDATA.W = OIDATA.W[:, OIDATA.indwvl]

    OIDATA.T3 = OIDATA.T3[:, OIDATA.indwvl] * unitconvert
    OIDATA.T3err = OIDATA.T3err[:, OIDATA.indwvl] * unitconvert
    if OIDATA.isDP == 1
        OIDATA.DP = OIDATA.DP[:, OIDATA.indwvl] * unitconvert
        OIDATA.DPerr = OIDATA.DPerr[:, OIDATA.indwvl] * unitconvert
    end
    OIDATA.wvl = vec(mean(OIDATA.wvl[:, OIDATA.indwvl], 1))
# # compute Differential Phases as needed for PAINTER and adjust angle unit
    OIDATA.nb = size(OIDATA.U, 1)
    OIDATA.nw = size(OIDATA.U, 2)

# # Create cluster of independent CLosure Phases and associated
    OIDATA = independentT3data(OIDATA)

    return OIDATA
end

function independentT3data(OIDATA::PAINTER_Input)
  Closure_index = OIDATA.Closure_index
  Cluster = independentT3(Closure_index)
  baseNb = basesincluster(Cluster)
  orderedCluster = makeclusterordered(Cluster,baseNb)
  rowt3 = t3row(Closure_index,Cluster)
  Nt3indep = length(baseNb)
  Xi = Dict{}()
  K = Dict{}()
  H = Dict{}()
  for n in 1:Nt3indep
      T3 = OIDATA.T3[rowt3[n],:]
      T3err = OIDATA.T3err[rowt3[n],:]
      if OIDATA.isDP == 1
          HDP = phasetophasediff(orderedCluster[n], OIDATA.nw, length(baseNb[n]), 0, 1, OIDATA.dptype, OIDATA.dpprm)
          DP = OIDATA.DP[baseNb[n],:]
          DPerr = OIDATA.DPerr[baseNb[n],:]
          Xi[n] = vcat(vec(T3),HDP*vec(DP)  )
          K[n] = vcat(vec(T3err),abs.(HDP)*vec(DPerr)  )
      else
          Xi[n] = vec(T3)
          K[n] = vec(T3err)
      end
  end
  OIDATA.Xi = Xi
  OIDATA.K = K
  OIDATA.baseNb = baseNb
  OIDATA.orderedCluster= orderedCluster
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
    U_3D = repmat(U_Coord, 1, length(wvl)) ./ repmat(wvl', length(U_Coord), 1)
    V_3D = repmat(V_Coord, 1, length(wvl)) ./ repmat(wvl', length(V_Coord), 1)
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
    utmp = round.(1000 * V2_U_Coord )
    U1 = round.(1000 * T3_U1_Coord)
    U2 = round.(1000 * T3_U2_Coord)
    U3 = round.(1000 * (T3_U1_Coord + T3_U2_Coord))
    for k in 1:length(U1)
        pos1 = LT3 + find(x->(x == U1[k]), utmp)
        pos2 = LT3 + find(x->(x == U2[k]), utmp)
        pos3 = LT3 + find(x->(x == U3[k]), utmp)
        Closure_index = vcat(Closure_index, hcat(pos1, pos2, pos3))
    end
    LT3 = maximum(Closure_index[:])
    LFil = LFil + length(U1)
    return Closure_index, LT3, LFil
end
# ---------------------------------------------------------------------------------
# PAINTER used Differential Visibilities, creation from Differential Visibilities
# ---------------------------------------------------------------------------------
function diffphi(DiffPhi::Array, HDP::SparseMatrixCSC)
# DiffPhi: Differential visibilities as given in OIFITS Data
# DiffPhi is matrix (Nbase*Nwvl)
# All channels reference are assumed to be equal
# defined as in equation 21 of PAINTER (1)
    DiffPhiAB = HDP * vec(DiffPhi)
end
# ---------------------------
# All channels are assumed to be independent
# variance of sum of random variables = sum of random variables variance
function diffphierr(DiffPhierr::Array, HDP::SparseMatrixCSC)
    DiffPhiABerr = abs(HDP) * vec(DiffPhierr)
end
# ---------------------------------------------------------------------------------
# INDEPENDENT Phases Closures from index -- to put in PAINTEROIFITS.JL
# ---------------------------------------------------------------------------------
function independentT3(Closure_index::Matrix)
    Cluster = Dict{}()
    Index = collect(1:size(Closure_index,1))
    m=0
    while sum(Index)>0
        m+=1
        idxinit = find(Index.>0)[1]
        Cluster[m] = reshape(Closure_index[Index[idxinit],:],1,3)
        Index[idxinit] = 0
        nn = 0
        while true
            nn+=1
            init = vec(Cluster[m][nn,:])
            Cluster[m], Index = searchfromclusterinlist(Cluster[m],init,Index,Closure_index)
            l = size(Cluster[m],1)
            if  nn==l
                break
            end

          end
      end
    return Cluster
end
function searchfromclusterinlist(Cluster::Matrix,init::Vector,Index::Vector,Closure_index::Matrix)
    for n in Index
      if n>0
        if( !isempty( find( init[1] .== Closure_index[n,:] ) )
          ||!isempty( find( init[2] .== Closure_index[n,:] ) )
          ||!isempty( find( init[3] .== Closure_index[n,:] ) ) )
            Cluster = vcat( Cluster , reshape(Closure_index[n,:],1,3) )
            Index[n] = 0
        end
      end
    end
    return Cluster, Index
end
# ---------------------------------------------------------------------------------
# Number of bases in each cluster - tools to create smaller T3 and DP Matrices
# ---------------------------------------------------------------------------------
function basesincluster(Cluster::Dict)
    Nc = length(Cluster)
    basenb = Dict{}()
    for n in 1:Nc
        vecclust = sort( vec( Cluster[n] ) )
        indpos = vcat(1, diff( vecclust ).>0)
        basenb[n] = round.(Int, vecclust[ (vecclust.*indpos).>0 ])
    end
    return basenb
end
# ---------------------------------------------------------------------------------
# Sort bases in Cluster from 1 to Nb base in Cluster
# ---------------------------------------------------------------------------------
function makeclusterordered(Cluster::Dict,baseNb::Dict)
  Nc = length(Cluster)
  theorderedcluster = Dict{}()
  for n in 1:Nc
      thecluster = Cluster[n]
      baseinthecluster = baseNb[n]
      numberofindependentbaseinthecluster = length(baseinthecluster)
      nt3 = size(thecluster,1)
      orderedbaseinthecluster = collect(1:numberofindependentbaseinthecluster)
      thevectorizedcluster = vec(thecluster)
      theorderedvectorizedcluster = zeros(thevectorizedcluster)
      for m in 1:numberofindependentbaseinthecluster
          cond = find( thevectorizedcluster .==  baseinthecluster[m])
          theorderedvectorizedcluster[cond] = orderedbaseinthecluster[m]
        end
      theorderedcluster[n] = reshape(theorderedvectorizedcluster,nt3,3)
  end
  return theorderedcluster
end
# ---------------------------------------------------------------------------------
# T3 rows related to H
# ---------------------------------------------------------------------------------
function t3row(Closure_index::Matrix,Cluster::Dict)
  Nc = length(Cluster)
  rowt3 = Dict{}()
  for n in 1:Nc
      thecluster = Cluster[n]
      lc = size(thecluster,1)
      rowt3[n] = zeros(lc)
      for m in 1:lc
          rowt3[n][m] = find( prod(Closure_index .== reshape(thecluster[m,:],1,3),2))[1]
      end
      rowt3[n] = round.(Int,rowt3[n])
  end
  return rowt3
end
