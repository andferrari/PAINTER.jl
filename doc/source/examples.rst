.. _examples-label:

Examples and demo
=================

Demo for impatients
-------------------

``PAINTER.jl`` contains a demo file ``painterdemo.jl``
with an OIFITS folder in the default installation folder.
To run the demo type:

.. code:: julia

  using PAINTER
  painterdemo()

``painterdemo()`` run a simation with data generated with ASPRO with AMBER configuration and a gray object.
The demo includes warm start, save and load of structures, a custom plot function (require PyPLot), ...

``painterdemo("gravity")`` run simulation with data from the beauty contest 2016 (http://www.opticalinterferometry.com/beauty2016).
Data will be downloaded in current folder and contains ``gravity`` simulation.
In this case PAINTER uses the phases of the complexe visibities and the closure phases for the phases estimation.
The demo includes save and load of structures, a custom plot function (require PyPLot), ...

``painterdemo("bc04")`` run simulation with data from the beauty contest 2004.
Data are monochromatic. The demo includes save and load of structures, a custom plot function (require PyPLot), ...

User parameters and single execution for ``painterdemo()``
----------------------------------------------------------

* The folowing parameters are set by the user:

  .. code:: julia

    path        = '../OifitsFolder'
    FOV         = 0.01
    indwvl      = 1:30
    nx          = 64
    eps1        = 1e-4
    eps2        = 1e-4
    rho_y       = 10
    alpha       = 1e4
    beta        = 1e5
    rho_spat    = 4
    rho_ps      = rho_spat
    lambda_spat = 1e-5
    rho_spec    = 1/2
    lambda_spec = 1e-5
    dptype      = "sliding" # type of differential phases
    aff         = true      # plot is enabled
    nbitermax   = 100

  ``PAINTER.jl`` will extract OIFITS informations from all files in the folder ``../OifitsFolder`` and will restrict the analysis to the first 29 wavelengths.

* The initial estimate is the default.  ADMM is enabled by default and will run the algorithm for 100 iterations.
* The support constraint is defined by a disk:

  .. code:: julia

    mask3D = mask(nx, int(nx/2 -3), choice="disk")

* Other parameters take the default values.

``PAINTER.jl`` is then executed:

.. code:: julia

  OIDATA, PDATA = painter(Folder=Folder, nbitermax=nbitermax, nx=nx, lambda_spat=lambda_spat=lambda_spat, lambda_spec=lambda_spec, rho_y= rho_y, rho_spat= rho_spat, rho_spec= rho_spec, rho_ps= rho_ps, alpha= alpha, beta=beta, eps1=eps1, eps2=eps2, FOV= FOV, indwvl=indwvl)

Algorithm warm start
--------------------

``PDATA`` contains all variables and array modified during iterations, including the Lagrange
multipliers. This allows a warm start of the ADMM algorithm. This is useful for example when
the iterations have been stopped by ``nbitermax`` but the algorithm has not yet converged.

In this example the user wants 1000 additional iterations with disabled plots:

.. code:: julia

  nbitermax += 1000
  aff = false
  OIDATA, PDATA_new = painter(OIDATA,PDATA, nbitermax, aff, PlotFct = Plotfunction)

``PDATA_new`` is used to store the new auxiliary variables.

Outer iterations mode
---------------------

It is possible to save the estimates (or other variables) at each iteration
using single iterations in a loop:

.. code:: julia

    for n = 1:10
      nbitermax += 1
      OIDATA, PDATA = painter(OIDATA, PDATA, nbitermax, aff)
      saveX[n] = PDATA.x
      saveW[n] = PDATA.w
    end

Note that this is a very time consuming process.

User defined plot function
--------------------------

It is possible to plot or to print some informations on available data during iterations.
If ``PyPlot.jl`` is installed, ``painter`` will execute each ``CountPlot`` iterations the function defined by the variable ``PlotFct``. This user defined function must respect the input arguments of ``painterplotfct``:

.. function:: Plotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)

For example, to plot at each iteration the sum over all wavelengths of an estimated polychromatic  object, projected on a support constraint:

.. code:: julia

	using PyPlot

	function Plotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)
		x = PDATA.x
		s = (PDATA.w.>0.0)
		im2show = squeeze(sum(x.*s,3),3)
		imshow(im2show)
	end

	OIDATA,PDATA = painter(..., PlotFct = Plotfunction)
