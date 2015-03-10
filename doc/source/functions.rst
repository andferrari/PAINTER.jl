Functions
=========

Main functions
~~~~~~~~~~~~~~

.. function:: painter(...)

  ``painter`` is defined with two methods:

  -  full parameters definition. This method is generally used to initialize the algorithm:

    .. code:: julia

      OIDATA,PDATA,OPTOPT = painter(Folder, nbitermax, nx, lambda_spat, lambda_spec, lambda_L1, epsilon, rho_y, rho_spat, rho_spec, rho_ps, alpha, Wvlt, beta, eps1, eps2, FOV, mask3D, xinit3D, indfile, indwvl, ls, scl, gat, grt, vt, memsize, mxvl, mxtr, stpmn, stpmx, aff, CountPlot, admm, paral)

  -  Specific structures. This method allows to restart the algorithm,
     for example if the number of iterations is not sufficient (see variable ``nbitermax+=100``).

    .. code:: julia

        OIDATA,PDATA,OPTOPT = painter(OIDATA,PDATA,OPTOPT, nbitermax, aff)

  ``painter`` returns 3 structures:

  .. code:: julia

    OIDATA,PDATA,OPTOPT = painter(...)

  where:

  - ``OIDATA``: contains all oifits information and user defined parameters
  - ``PDATA``: contains all variables and array modified during iterations
  - ``OPTOPT``: contains all OptimPack parameters for the phases minimization process


.. function:: mask(nx::Int,param::Int,choice::ASCIIString)

  Creates a binary mask of size nx:sup:`2`:
E = mc\ :sup:`2`
  .. code:: julia Mymask3D = mask(nx,param,choice)

  - ``choice`` can be a square (default: ``choice="square"``) or a
  disk (``choice="disk"``).
  - ``nx`` is the size of the image.
  - ``param`` is the radius of the disk or the half size of the square.

.. function:: paintersave(savepath::ASCIIString,PDATA::PAINTER_Data,OIDATA::PAINTER_Input,OPTOPT::OptOptions)


  Save structures ``OIDATA``, ``PADATA`` and ``OPTOPT`` (TBD) into ``*.jld`` (see `HDF5<https://github.com/timholy/HDF5.jl>`_ package)
  files. The names of the files are defined by ``savepath``. To save the structures in ``file.jld``

  .. code:: julia

    savepath = "../Mypath/file.jld"
    paintersave(savepath,PDATA,OIDATA,OPTOPT)


.. function:: painterload(savepath::ASCIIString)


  To load the structures contained in ``*.jld`` files:

  .. code:: julia

    PDATA2,OIDATA2,OPTOPT2 = painterload(savepath)

painterplotfct function (painterplot.jl):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: julia

    function painterplotfct(x::SharedArray, w::Array, crit1::Vector, crit2::Vector, eps1::Real, eps2::Real, nx::Int64, nw::Int64, wvl::Vector, FOV::Real)

In order to allow user to draw personalized plots ``painterplot.jl`` is
a separated files of the package located in
``Painter.jl/src/painterplot.jl``. The default function compute
automatically number of subplot as a function of the number of
wavelength (if nw<30) and draw on the first figure the per-channel
estimates projected on the positiv support constraint. The axis are
defined by the field of view with no limitation of the amplitude
(colorbars are different for all images). A second figure draw the
primal and dual residuals (``crit1``\ and ``crit2``) as a function of
the iteration number and print the verbose of these values.
