# Structure for VMLM of OptimPack, see OptimPack.jl
# https://github.com/emmt/OptimPack.jl

ls = OptimPack.MoreThuenteLineSearch(ftol = 1e-8, gtol = 0.95)
scl = OptimPack.SCALING_OREN_SPEDICATO
gat = 0
grt = 1e-3
vt = false
memsize = 100
mxvl = 1000
mxtr = 1000
stpmn = 1e-20
stpmx = 1e+20
