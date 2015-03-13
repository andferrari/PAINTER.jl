.. _examples-label:

Examples
========


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


Demo code for impatients
------------------------

The following code can be used to test all the functionalities of the algorithm. This demo will reconstruct a nx\ :sup:`2` pixels gray object. The data are stored in 4 OIFITS files. They result from 102 bases at 227 wavelength and 34 phases closure per wavelength. As on the first example, the analysis is done on the first 29 wavelengths using all files, the field of view is 0.01 arc second. The execution will be parallelized and at each 10 iterations the 29 estimates of the object will be plot is ``PyPlot`` is installed.

.. code:: julia


	using PyPlot

	function myplotfunction(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)
		nx = OIDATA.nx
		nw = OIDATA.nw
		wvl = OIDATA.wvl
		FOV = OIDATA.FOV
		x = PDATA.x
		w = PDATA.w .> 0.

		indpix = linspace(-(FOV / 2), (FOV / 2), nx)
		pos = int([1, round(nx / 4), round(nx / 2), round(nx * 3 / 4), nx])

		count_y = 0
		count_x = 0
		SubRow  = 6
		SubColumn = 5

		for n in 1:nw
			subplot(SubColumn, SubRow, n)
			imshow(x[:, :, n] .* max(0, w[:, :, n]), origin = "lower")
			titlestring = @sprintf("%2.4f µm", wvl[n] * 1e6)
			title(titlestring)
			xticks([])
			yticks([])
			if( n == (nw + 1 - SubRow + count_x) )
				xticks([pos - 1], round(indpix[pos] * 100000) / 100)
				xlabel("FOV (mas)")
				count_x += 1
			end
			if(n == (1 + count_y * SubRow))
				yticks([pos - 1], round(indpix[pos] * 100000) / 100)
				ylabel("FOV (mas)")
				count_y += 1
			end
		end
	end

    MyPlotFct = myplotfunction
    MyFOV = 0.01
    Myindwvl = 1:30
    Mynx = 64
    Myeps1 = 1e-4
    Myeps2 = 1e-4
    Myrho_y = 10
    Myalpha = 1e4
    Mybeta = 1e5
    Myrho_spat = 4
    Myrho_ps = Myrho_spat
    Mylambda_spat = 1e-5
    Myrho_spec = .5
    Mylambda_spec = 1e-5
    Myaff = true     # plot is enabled
    Mynbitermax = 1
    Mypar = true     # parallel computing is disabled
    savepath = "mydata.jld"

    # initialize algorithm and run admm
    OIDATA, PDATA, OPTOPT = painter(nbitermax = Mynbitermax, nx = Mynx, lambda_spat = Mylambda_spat,
                                lambda_spec = Mylambda_spec, rho_y = Myrho_y, rho_spat = Myrho_spat,
                                rho_spec = Myrho_spec, rho_ps = Myrho_ps, alpha = Myalpha, beta = Mybeta,
                                eps1 = Myeps1, eps2 = Myeps2, FOV = MyFOV, indwvl = Myindwvl,
                                ls = OptimPack.MoreThuenteLineSearch(ftol = 1e-8, gtol = 0.95),
                                scl = OptimPack.SCALING_OREN_SPEDICATO, gat = 0, grt = 1e-3,
                                vt = false, memsize = 100, mxvl = 1000, mxtr = 1000, stpmn = 1e-20,
                                stpmx = 1e+20, PlotFct = MyPlotFct, aff = Myaff)


    # save data struture in .jld files
    paintersave(savepath,PDATA,OIDATA,OPTOPT)

    # load data struture in .jld files
    PDATA, OIDATA = painterload(savepath)

    # Warm start of the algorithm
    OIDATA, PDATA, OPTOPT = painter(PDATA,OIDATA,OPTOPT,100,true, PlotFct = MyPlotFct)
