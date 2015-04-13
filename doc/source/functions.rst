Functions
=========

Main function
-------------

.. function:: painter(...)

* ``painter`` is defined with two methods:

  * Full parameters definition. This method is generally used to initialize the algorithm:

    .. code:: julia

      OIDATA,PDATA,OPTOPT = painter(Folder, nbitermax, nx, lambda_spat, lambda_spec, lambda_L1, epsilon, rho_y, rho_spat, rho_spec, rho_ps, alpha, Wvlt, beta, eps1, eps2, FOV, mask3D, xinit3D, indfile, indwvl, ls, scl, gat, grt, vt, memsize, mxvl, mxtr, stpmn, stpmx, aff, CountPlot,PlotFct, admm, paral)

  * Specific structures. This method allows to restart the algorithm, for example if the number of iterations is not sufficient (see variable ``nbitermax+=100``).

    .. code:: julia

        OIDATA,PDATA,OPTOPT = painter(OIDATA,PDATA,OPTOPT, nbitermax, aff; plotfunction)

* ``painter`` returns 3 structures:

  .. code:: julia

    OIDATA,PDATA,OPTOPT = painter(...)

  where:

  * ``OIDATA``: contains all OIFITS information and user defined parameters.
  * ``PDATA``: contains all the variables and arrays modified during iterations.
  * ``OPTOPT``: contains all OptimPack parameters for the phases proximal operator.

Auxiliary functions
-------------------

.. function:: paintersave(savepath::ASCIIString,PDATA::PAINTER_Data,OIDATA::PAINTER_Input,OPTOPT::OptOptions)

  Saves the structures ``OIDATA``, ``PDATA`` and ``OPTOPT`` (TBD) into ``*.jld`` julia data files. See `HDF5 <https://github.com/timholy/HDF5.jl>`_ package.

  .. code:: julia

    savepath = "../Mypath/file.jld"
    paintersave(savepath,PDATA,OIDATA,OPTOPT)

.. function:: painterload(savepath::ASCIIString)


  Loads the structures contained in ``*.jld`` files:

  .. code:: julia

    PDATA2,OIDATA2 = painterload(savepath)

  The current version of the save function cannot save, and so load, the Optim structure OPTOPT and the pointer to the user defined plot function. To warmstart the algorithm, the user must redefine an OPTOPT structure and call the ``painter(...)`` with the personalized plot function as argument otherwise the default plot function is used.	

.. function:: painterfitsexport(savepath::ASCIIString,PDATA::PAINTER_Data, OIDATA::PAINTER_Input; forceWvlExt=false)

  Saves the relevant information from  ``PDATA`` (ouput data cube and associated criteria) and from  ``OIDATA`` (wavelengths, input reconstruction parameters,...) into a FITS file. The wavelengths are either implicitely defined in the main header (as for x and y axes) if they have a regular spacing, or explicitely given in a dedicated binary table, if they are irregular or if keyword parameter forceWvlExt=true.

  .. code:: julia

    savepath = "../Mypath/myfitsFile.fits"
    painterfitsexport(savepath,PDATA,OIDATA)

	
.. function:: painterplotfct(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)

  It is recommended to monitor the iterations of the algorithm when the number
  of wavelength is small, e.g. < 30.

  The default function computes the number of subplots as a function of the number of wavelength if ``nw<30``.
  Its is automatically called if ``PyPlot`` is installed and ``aff=true``.

* The first figure shows the per-channel estimates projected on the domain support. The axis are defined by the field of view with no limitation of the amplitude (colorbars are different for all images).
* The second figure shows the primal and dual residuals (``crit1`` and ``crit2``) as a function of the iteration.

  
.. function:: mask(nx::Int,param::Int,choice::ASCIIString)

  Creates a binary mask of size nx\ :sup:`2`:

    .. code:: julia

      Mymask3D = mask(nx,param,choice)

* ``choice`` can be a square (default: ``choice="square"``) or a disk (``choice="disk"``).
* ``nx`` is the size of the image.
* ``param`` is the radius of the disk or the half size of the square.
