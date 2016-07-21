# using PAINTER
using Base.Test

FOV = 0.01
indwvl = 1:30
nx = 64
eps1 = 1e-4
eps2 = 1e-4
rho_y = 10.
alpha = 1e4
beta = 1e5
rho_spat = 4.
rho_ps = rho_spat
lambda_spat = 1e-5
rho_spec = 0.5
lambda_spec = 1e-5
nbitermax = 4
dptype = "sliding"
dpprm = 5
aff=false

OIDATA, PDATA = painter(nbitermax = nbitermax, nx = nx, lambda_spat = lambda_spat,
                            lambda_spec = lambda_spec, rho_y = rho_y, rho_spat = rho_spat,
                            rho_spec = rho_spec, rho_ps = rho_ps, alpha = alpha, beta = beta,
                            eps1 = eps1, eps2 = eps2, FOV = FOV, indwvl = indwvl,
                            dptype = dptype, dpprm = dpprm, aff=aff,pathoptpkpt="")


@test_approx_eq_eps PDATA.crit1[1] 1.92717 1E-4
