.. _examples-label:

Examples and demo
=================

Demo for impatients
-------------------

``Painter.jl`` contains a demo file ``painterdemo.jl``
with an OIFITS folder in the default installation folder.
To run the demo type:

.. code:: julia

  demo = string(Pkg.dir("Painter"),"/src/painterdemo.jl")
  include(demo);

The demo includes warm start, save and load of structures, a custom plot function, ...

User parameters and single execution
------------------------------------

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

  ``Painter.jl`` will extract OIFITS informations from all files in the folder ``../MyOifitsFolder`` and will restrict the analysis to the first 29 wavelengths.

* The initial estimate is the default.  ADMM is enabled by default and will run the algorithm for 100 iterations.
* The support constraint is defined by a disk:

  .. code:: julia

    Mymask3D = mask(Mynx, int(Mynx/2 -3), "disk")

* Other parameters take the default values.

``Painter.jl`` is then executed:

.. code:: julia

  OIDATA, PDATA, OPTOPT = painter(Folder=MyFolder, nbitermax=Mynbitermax, nx=Mynx, lambda_spat=Mylambda_spat=Mylambda_spat, lambda_spec=Mylambda_spec, rho_y= Myrho_y, rho_spat= Myrho_spat, rho_spec= Myrho_spec, rho_ps= Myrho_ps, alpha= Myalpha, beta=Mybeta, eps1=Myeps1, eps2=Myeps2, FOV= MyFOV, indwvl=Myindwvl, paral=Myparal)

Algorithm warm start
--------------------

``PDATA`` contains all variables and array modified during iterations, including the Lagrange
multipliers. This allows a warm start of the ADMM algorithm. This is usefull for example when
the iterations have been stopped by ``nbitermax`` but the algorithm has not yet converged.

In this example the user wants 1000 aditional iterations with disabled plots:

.. code:: julia

  nbitermax += 1000
  aff = false
  OIDATA, PDATA_new, OPTOPT = painter(OIDATA,PDATA,OPTOPT, nbitermax, aff, PlotFct = myPlotfunction)

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

User defined plot function
--------------------------

It is possible to plot or to print some informations on available data during iterations.
If ``PyPlot.jl`` is installed, ``painter`` will execute each ``CountPlot`` iterations the function defined by the variable ``PlotFct``. This user defined function must respect the input arguments of ``painterplotfct``:

.. function:: myPlotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)

For example, to plot at each iteration the sum over all wavelengths of an estimated polychromatic  object, projected on a support constraint:

.. code:: julia

	using PyPlot

	function myPlotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)
		x = PDATA.x
		s = (PDATA.w.>0.0)
		im2show = squeeze(sum(x.*s,3),3)
		imshow(im2show)
	end

	OIDATA,PDATA,OPTOPT = painter(..., PlotFct = myPlotfunction)
