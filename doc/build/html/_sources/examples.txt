Examples
========


User parameters and single algorithm execution
----------------------------------------------

* The folowing parameters are set by the user:

  .. code:: julia

    Mypath        = '../MyOifitsFolder'
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
    Myaff         = true      # plot is enabled
    Mynbitermax   = 100
    Mypar         = false     # parallel computing is disabled

  ``Painter.jl`` will extract OIFITS informations from all files in ``Mypath`` and will restrict the analysis to the first 29 wavelengths.

* The initial estimate is the default.   admm is enabled by default and will run the algorithm for 100 iterations.
* The support constraint is defined by a disk:

  .. code:: julia

    Mymask3D = mask(Mynx, int(Mynx/2 -3), "disk")

* Other parameters take the default values.

* ``Painter.jl`` is then executed:

  .. code:: julia

    OIDATA, PDATA, OPTOPT = painter(Folder=MyFolder, nbitermax=Mynbitermax, nx=Mynx, lambda_spat=Mylambda_spat=Mylambda_spat, lambda_spec=Mylambda_spec, rho_y= rho_y, rho_spat= rho_spat, rho_spec= rho_spec, rho_ps= rho_ps, alpha= alpha, beta=Mybeta, eps1=Myeps1, eps2=Myeps2, FOV= MyFOV, indfile, indwvl=Myindwvl, paral=Myparal)

Algorithm warmstart
-------------------

``PDATA`` contains all variables and array modified during iterations, including the Lagrange
multipliers. This allows a warmstart of the ADMM algorithm. This is usefull for example when
the iterations have been stopped by ``nbitermax`` but the algorithm has not converged.

In this example the user wants 1000 aditional iterations with disabled plots:

.. code:: julia

  nbitermax += 1000
  aff = false
  OIDATA, PDATA_new, OPTOPT = painter(OIDATA,PDATA,OPTOPT, nbitermax, aff)

``PDATA_new`` is used to store the new auxiliary variables.

Outer iterations mode
---------------------

It is possible to save the estimates (or other variables) at each iteration
using single iterations in a loop:

.. code:: julia

    for n = 1:10
      nbitermax += 1
      OIDATA, PDATA, OPTOPT = painter(OIDATA, PDATA, OPTOPT, nbitermax, aff)
      saveX[n] = PDATA.x
      saveW[n] = PDATA.w
    end

Note that this is a very time consuming process.
