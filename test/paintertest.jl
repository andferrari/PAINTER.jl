using PAINTER
using Base.Test

MyFOV = 0.01
Myindwvl = 1:30
Mynx = 64
Myeps1 = 1e-4
Myeps2 = 1e-4
Myrho_y = 10.
Myalpha = 1e4
Mybeta = 1e5
Myrho_spat = 4.
Myrho_ps = Myrho_spat
Mylambda_spat = 1e-5
Myrho_spec = 0.5
Mylambda_spec = 1e-5
Mynbitermax = 4

OIDATA, PDATA, OPTOPT = painter(nbitermax = Mynbitermax, nx = Mynx, lambda_spat = Mylambda_spat,
                            lambda_spec = Mylambda_spec, rho_y = Myrho_y, rho_spat = Myrho_spat,
                            rho_spec = Myrho_spec, rho_ps = Myrho_ps, alpha = Myalpha, beta = Mybeta,
                            eps1 = Myeps1, eps2 = Myeps2, FOV = MyFOV, indwvl = Myindwvl,
                            ls = OptimPack.MoreThuenteLineSearch(ftol = 1e-8, gtol = 0.95),
                            scl = OptimPack.SCALING_OREN_SPEDICATO, gat = 0, grt = 1e-3,
                            vt = false, memsize = 100, mxvl = 1000, mxtr = 1000, stpmn = 1e-20,
                            stpmx = 1e+20);


@test_approx_eq_eps PDATA.crit1[1] 1.92813 1E-4
