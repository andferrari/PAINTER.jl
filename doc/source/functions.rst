Functions
=========

Main function
-------------

.. function:: painter(...)

* ``painter`` is defined with two methods:

  * Full parameters definition. This method is generally used to initialize the algorithm:

    .. code:: julia

      OIDATA,PDATA = painter(Folder, nbitermax, nx, lambda_spat, lambda_spec, lambda_L1, epsilon, rho_y, rho_spat, rho_spec, rho_ps, alpha, Wvlt, beta, eps1, eps2, FOV, mask3D, xinit3D, indfile, indwvl, aff, CountPlot, PlotFct, admm)

  * Specific structures. This method allows to restart the algorithm, for example if the number of iterations is not sufficient (see variable ``nbitermax+=100``).

    .. code:: julia

        OIDATA,PDATA = painter(OIDATA,PDATA, nbitermax, aff; plotfunction)

* ``painter`` returns 3 structures:

  .. code:: julia

    OIDATA,PDATA = painter(...)

  where:

  * ``OIDATA``: contains all OIFITS information and user defined parameters.
  * ``PDATA``: contains all the variables and arrays modified during iterations.

Auxiliary functions
-------------------

.. function:: paintersave(filename::String,PDATA::PAINTER_Data,OIDATA::PAINTER_Input)

  Saves the structures ``OIDATA``, ``PDATA`` into ``*.jld`` julia data files. The prefix of these structures is added before the "filename" base when writing the output files. See `HDF5 <https://github.com/timholy/HDF5.jl>`_ package for details on the format.

  .. code:: julia

    filename= "datafile.jld"
    cd("~/path/to/saved/data") # move to a different directory if necessary
    paintersave(filename,PDATA,OIDATA)

.. function:: painterload(filename::String)


  Loads the structures from ``*.jld`` files. The files to be loaded must start with OIDATA_ and PDATA_ prefixes, but the filename entered as an argument should not have a prefix, since they are internally added by this function. Therefore, the filename of ``painterload`` is compatible with the one of ``paintersave``.

  .. code:: julia

    PDATA2,OIDATA2 = painterload(filename)

  The current version of the save function doesn't save the pointer to the user defined plot function. To warmstart the algorithm, the user must call the ``painter(...)`` with the personalized plot function as argument otherwise the default plot function is used.

.. function:: painterfitsexport(filename::String,PDATA::PAINTER_Data, OIDATA::PAINTER_Input)

  Saves the relevant information from  ``PDATA`` (output data cube and associated criteria, reconstructed complex visibilities,...) and from  ``OIDATA`` (wavelengths, input reconstruction parameters,...) into a FITS file "filename", which possibly includes a full path. The resulting FITS file has three HDUs : "Primary" is the reconstructed image cube, "INFO" contains the reconstruction parameters and criteria, and "VIS" contains the complex visibilities of the reconstruction, with the associated wavelengths and (U,V) points.

  .. code:: julia

    filename = "~/path/to/saved/data/myfitsdata.fits"
    painterfitsexport(filename,PDATA,OIDATA)


.. function:: painterplotfct(PDATA::PAINTER_Data,OIDATA::PAINTER_Input)

  It is recommended to monitor the iterations of the algorithm when the number
  of wavelength is small, e.g. < 30.

  The default function computes the number of subplots as a function of the number of wavelength if ``nw<30``.
  Its is automatically called if ``PyPlot`` is installed and ``aff=true``.

* The first figure shows the per-channel estimates projected on the domain support. The axis are defined by the field of view with no limitation of the amplitude (colorbars are different for all images).
* The second figure shows the primal and dual residuals (``crit1`` and ``crit2``) as a function of the iteration.


.. function:: mask(nx::Int,param::Int,choice::String)

  Creates a binary mask of size nx\ :sup:`2`:

    .. code:: julia

      Mymask3D = mask(nx,param,choice)

* ``choice`` can be a square (default: ``choice="square"``) or a disk (``choice="disk"``).
* ``nx`` is the size of the image.
* ``param`` is the radius of the disk or the half size of the square.
