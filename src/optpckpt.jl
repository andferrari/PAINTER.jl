# Structure for VMLM of OptimPack, see OptimPack.jl
# https://github.com/emmt/OptimPack.jl
ls = OptimPack.MoreThuenteLineSearch(ftol = 1e-8, gtol = 0.95)
gat = 0
grt = 1e-6
vt = false
memsize = 10
mxvl = 1000
mxtr = 1000
