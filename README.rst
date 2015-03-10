PAINTER.jl
==========

|Build Status|

``PAINTER.jl`` is a julia implementation of PAINTER: Polychromatic
opticAl INTErferometric Reconstruction software described in [1] and
[2].

Installation
------------

``PAINTER.jl`` uses the library
```OptimPack`` <https://github.com/emmt/OptimPack>`__ for solving the
phases proximal operator. The Julia interface to OptimPack is not a
registered Julia package. To install
```OptimPack.jl`` <https://github.com/emmt/OptimPack.jl>`__ type from a
Julia session the following commands:

.. code:: julia

    Pkg.clone("https://github.com/emmt/OptimPack.jl.git")
    Pkg.build("OptimPack")

Installation of PAINTER is then as simple as typing:

.. code:: julia

    Pkg.clone("https://github.com/andferrari/PAINTER.jl.git")

Dependencies
------------

PAINTER uses the following registered Julia package. They will be
automaticaly installed during PAINTER installation.

-  ```FITSIO.jl`` <https://github.com/JuliaAstro/FITSIO.jl>`__: Julia
   support for OI-FITS (optical interferometry data format).
-  ```NFFT.jl`` <https://github.com/tknopp/NFFT.jl>`__: Julia
   implementation of the Non-equidistant Fast Fourier Transform (NFFT).
-  ```Wavelets`` <https://github.com/JuliaDSP/Wavelets.jl>`__: A Julia
   package for fast wavelet transforms.
-  ```HDF5`` <https://github.com/timholy/HDF5.jl>`__: for writing JLD
   ("Julia data") variables.

It is recommended to install
```PyPlot.jl`` <https://github.com/stevengj/PyPlot.jl>`__ to monitor the
iterations of the algorithm. See
```PyPlot.jl`` <https://github.com/stevengj/PyPlot.jl>`__ page.

Typical usage
-------------

``PAINTER.jl`` can be used with all parameters defined by function ###
painter function

``painter``\ return 3 structures:

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

Two parameters are independent of structures:

-  ``nbitermax``: number of
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__
   iteration, default 1000
-  ``aff``: true for plot drawing using ``PyPlot.jl``, default: false

where the folowing arguments, contained in ``OIDATA`` structure, define:

-  ``Folder``: path to the folder containing oifits/fits files (need
   ``OIFITS.jl``), default: src/OIFITS
-  ``nx`` : size in pixels of image (image of size nx\*nx), default 64
   64
-  ``Wvlt`` : list of wavelets basis, default: 8 Daubechies and haar
   wavelets
-  ``lambda_spat`` : Spatial Regularization parameter (weight) (Eq.29,31
   [1]), default: nx^-2
-  ``lambda_spec`` : Spectral Regularization parameter (weight)
   (Eq.29,31 [1]), default: 1e-2
-  ``lambda_L1`` : L1 constraint, emphasize sparsity of object (as stars
   field), default: 0
-  ``epsilon`` : Ridge/tikhonov parameter epilson \|X\|^2\_2 (Eq.29,31
   [1]), default 1e-6
-  ``rho_y`` :
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__
   parameter for data fidelity (convergence rate) (Eq.35,50-52 [1]),
   default: 1
-  ``rho_spat`` :
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__
   parameter for Spatial Regularization (convergence rate) (Eq.25,31
   [1]), default: 1
-  ``rho_spec`` :
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__
   parameter for Spectral Regularization (convergence rate) (Eq.42,55
   [1]), default: 1
-  ``rho_ps`` :
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__
   parameter for positivity (convergence rate) (Eq.47,54 [1]), default:
   1
-  ``alpha`` : Weight for Complexe visibility estimation from V2
   (Eq.25,31 [1]), default: 1
-  ``beta`` : Weight for Complexe visibility estimation from phases
   difference (Eq.25,31 [1]), default: 1
-  ``eps1`` : Primal Residual stop criterium (3.3 Optimality Conditions
   and Stopping Criterion , see
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__),
   default: 1e-6
-  ``eps2`` : Dual Residual stop criterium (3.3 Optimality Conditions
   and Stopping Criterion, see
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__),
   default: 1e-6
-  ``FOV`` : Field Of View User parameter, must be in ArcSecond, ,
   default: 40 mas
-  ``mask3D``: Support constraint, can be a path to a fits file, an user
   Array, an empty Array is default and returns a unconstraint support,
   or can be initialized with the ``mask(...)`` function. In case of
   array dimensions are checked and corrected, default: no constraint
-  ``xinit3D``: Initial Estimate, as for the mask, this parameter accept
   a path to a fits files. Can also be an Array, or complexe
   visibilities. In case of array dimensions are checked and corrected,
   default: centered dirac
-  ``indfile``: allow to reduce the set of file present in the
   ``Folder``, , default: all files
-  ``indwvl``: allow to restrict the set of wavelength used (wavelength
   index), default: all wavelengths
-  ``CountPlot``: draw plot at each ``CountPlot`` iterations, default:
   10
-  ``admm``: to run admm or not (true or false), false just initialize
   and give structures to run ``painteradmm(...)``, default: true
-  ``paral``: if true perform some calculus in parallel, useful when
   julia use several core (``julia -p nprocs``), default: true

``OPTOPT`` structure contains:

-  

   -  ``ls,scl,gat,grt,vt,memsize,mxvl,mxtr,stpmn,stpmx``: related to
      ```OptimPack`` <https://github.com/emmt/OptimPack>`__. default:
      ls=OptimPack.MoreThuenteLineSearch(ftol=1e-4,gtol=0.9),
      scl=OptimPack.SCALING\_OREN\_SPEDICATO, gat=1E-6, grt=1E-6,
      vt=false, memsize=100, mxvl=1000, mxtr=1000, stpmn=1E-20,
      stpmx=1E+20

if parameters are not setted, default value are used. For example,
calling: ``OIDATA,PDATA,OPTOPT =  painter()`` execute the 3D image
reconstruction algorithm from data stored in all \*.oifits files from
folder "OIFITS" located in the ``PAINTER.jl`` source folder
(``src/OIFITS``). The parameters are setted to default value with no
support contraint, spatial and spectral regularizations, positivity
constraint, the original estimate is a centered dirac at all
wavelengths.

Useful Array in ``PDATA`` are

-  ``PDATA.x``: reconstruced 3D images
-  ``PDATA.w``: positivity+support contraint, ``PDATA.x.*(PDATA.w.>0)``
   will project the reconstructed 3D images on positif support
-  ``PDATA.Fx``: non uniform Fourier transform of 3D images
-  ``PDATA.crit1`` and ``PDATA.crit2`` the primal and dual residual
   values of the
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__
   algorithm
-  ``PDATA.ind`` the iteration indice, useful to re-run algorithm

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

painterplot.jl: painterplotfct function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Examples:
~~~~~~~~~

User parameters and single algorithm execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The folowing parameters are setted byt the user, the initial estimate is
the default, drawing is enabled and parallel computing is disabled.
Painter will take oifits informations from all files and will restrict
the analysis on the first 29 wavelength. admm is enabled by default and
will run the algorithm for 1000 iterations.

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
    Myrho_ss      = 1
    Mylambda_spec = 1e-5
    Myaff         = true
    Mynbitermax   = 1000                 
    Mypar         = false

The support constraint is defined by a circle

-  ``Mymask3D = mask(Mynx,int(Mynx/2 -3),"circ")``

The other parameters are setted with default values. Painter is then
executed:

.. code:: julia

    OIDATA,PDATA,OPTOPT = painter(Folder=MyFolder, nbitermax, nx=Mynx, lambda_spat=Mylambda_spat=Mylambda_spat, lambda_spec=Mylambda_spec, rho_y= rho_y, rho_spat= rho_spat, rho_spec= rho_spec, rho_ps= rho_ps, alpha= alpha, beta=Mybeta, eps1=Myeps1, eps2=Myeps2, FOV= MyFOV, indfile, indwvl=Myindwvl, paral=Myparal)  

Algorithm re-run
^^^^^^^^^^^^^^^^

One can want to re-run the algorithm because the maximum number of
iteration is already reached and was not enough. In this example
consider that the algorithm was on the good way to converge to a good
solution, so to decrease the time cost of the algorithm the drawing is
disabled and to do another more 1000 iterations more, just set:

-  ``nbitermax += 1000``
-  ``aff="false"``

In order to keep last results another ouput structure is used
(``PDATA_new``)

Then re-run the algorithm:

.. code:: julia

    OIDATA,PDATA_new,OPTOPT = painter(OIDATA,PDATA,OPTOPT, nbitermax, aff)  

Iteration per Iteration reconstruction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the user desire to keep all the estimate (time consuming process),
iteration per iteration, consider to make a loop:

::

    for n=1:10
    nbitermax+=1
    OIDATA,PDATA,OPTOPT = painter(OIDATA,PDATA,OPTOPT, nbitermax, aff)  
    saveX[n] = PDATA.x
    saveW[n] = PDATA.w
    end

Authors
-------

PAINTER was developped at Laboratoire J.-L. Lagrange, Université de Nice
Sophia, CNRS, Observatoire de la Côte d'Azur, by `Antony
Schutz <http://www.antonyschutz.com>`__ and `André
Ferrari <https://www-n.oca.eu/aferrari>`__.

References
----------

The PAINTER algorithm is described in [1]. The original MATLAB code is
available `here <https://www-n.oca.eu/aferrari/painter/>`__ but the use
of ```PAINTER.jl`` <https://github.com/andferrari/PAINTER.jl>`__ is
highly recommended.
```PAINTER.jl`` <https://github.com/andferrari/PAINTER.jl>`__ implements
an accelerated version of PAINTER described in [2].

1. Schutz, A., Ferrari, A., Mary, D. Soulez, F., Thiébaut, E., Vannier,
   M. “PAINTER: a spatio-spectral image reconstruction algorithm for
   optical interferometry”. JOSA A. Vol. 31, Iss. 11, pp. 2356–2361
   (2014). `arxiv <http://arxiv.org/abs/1407.1885>`__
2. Schutz, A., Ferrari, A., Mary, D., Thiébaut, E., Soulez, F. “Large
   scale 3D image reconstruction in optical interferometry”. Submitted
   to EUSIPCO 2015, Nice. `arxiv <http://arxiv.org/abs/1503.01565>`__

Credits
-------

The development of OptimPack was partially supported by the
`POLCA <http://polca.univ-lyon1.fr>`__ project funded by the French
Agence Nationale pour la Recherche (ref. ANR-10-BLAN-0511).

License
-------

PAINTER is released under under the MIT "Expat" License.

.. |Build Status| image:: https://travis-ci.org/andferrari/PAINTER.jl.svg?branch=master
   :target: https://travis-ci.org/andferrari/PAINTER.jl
