functions
=========

painter functions
~~~~~~~~~~~~~~~~~

``PAINTER.jl`` can be used with all parameters defined by ``painter function``
and return 3 structures:

``OIDATA,PDATA,OPTOPT = painter(...)``

where:

-  ``OIDATA``: contains all oifits information and user defined
   parameters
-  ``PDATA``: contains all variables and array modified during
   iterations
-  ``OPTOPT``: contains all OptimPack parameters for the phases
   minimization process

it follows that ``painter`` can be used with two methods:

-  full parameters definition for first run algorithm:

.. code:: julia

    OIDATA,PDATA,OPTOPT = painter(Folder, nbitermax, nx, lambda_spat, lambda_spec, lambda_L1, epsilon, rho_y, rho_spat, rho_spec, rho_ps, alpha, Wvlt, beta, eps1, eps2, FOV, mask3D, xinit3D, indfile, indwvl, ls, scl, gat, grt, vt, memsize, mxvl, mxtr, stpmn, stpmx, aff, CountPlot, admm, paral)  

-  Given structures to re-run algorithm, for example if the defined
   number of iteration is not enough (nbitermax+=100)


.. code:: julia

    OIDATA,PDATA,OPTOPT = painter(OIDATA,PDATA,OPTOPT, nbitermax, aff)  

mask function
~~~~~~~~~~~~~

For an image of size nx^2, the support constraint binary mask can be
generated using the function:

``Mymask3D = mask(nx,param,choice)``

choice can be rectangular constraint (default: choice="rect") or a
circle (choice="circ"). nx is the side size of the image to reconstruct.
param is the radius of the circle or the half side of the square.

save and load
~~~~~~~~~~~~~

Thanks to ```HDF5`` <https://github.com/timholy/HDF5.jl>`__ package,
``save`` and ``load`` functions are provided with PAINTER.jl. To save
structures ``OIDATA``, ``PADATA`` and ``OPTOPT`` (TBD) into ``*.jld``
files defined by ``savepath``, its full path consists to call the
``save``\ function:

.. code:: julia

    function paintersave(savepath::ASCIIString,PDATA::PAINTER_Data,OIDATA::PAINTER_Input,OPTOPT::OptOptions)

example:
^^^^^^^^

to save the structures in ``file.jld``

.. code:: julia

    savepath = "../Mypath/file.jld"
    paintersave(savepath,PDATA,OIDATA,OPTOPT)

to load the structures into new structures

.. code:: julia

    PDATA2,OIDATA2,OPTOPT2 = painterload(savepath)

painterplotfct function (painterplot.jl):  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: julia

    function painterplotfct(x::SharedArray, w::Array, crit1::Vector, crit2::Vector, eps1::Real, eps2::Real, nx::Int64, nw::Int64, wvl::Vector, FOV::Real)  

In order to allow user to draw personalized plots ``painterplot.jl`` is
a separated files of the package located in
``PAINTER.jl/src/painterplot.jl``. The default function compute
automatically number of subplot as a function of the number of
wavelength (if nw<30) and draw on the first figure the per-channel
estimates projected on the positiv support constraint. The axis are
defined by the field of view with no limitation of the amplitude
(colorbars are different for all images). A second figure draw the
primal and dual residuals (``crit1``\ and ``crit2``) as a function of
the iteration number and print the verbose of these values.