variables
=========

``Painter.jl`` can be used with all parameters defined by ``painter function``
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

Variables and structures informations
-------------------------------------

Two parameters are independent of structures:

-  ``nbitermax``: number of
   ```admm`` <http://stanford.edu/~boyd/papers/admm_distr_stats.html>`__
   iteration, default 1000
-  ``aff``: true for plot drawing using ``PyPlot.jl``, default: false

User and data Structures
------------------------

where the folowing arguments, contained in ``OIDATA`` structure, define:

-  ``Folder``: path to the folder containing oifits/fits files (need
   ``OIFITS.jl``), default: src/OIFITS
-  ``nx`` : size in pixels of image (image of size nx\*nx), default: 64
-  ``Wvlt`` : list of wavelets basis, default: 8 Daubechies and haar
   wavelets
-  ``lambda_spat`` : Spatial Regularization parameter (weight) (Eq.29,31
   [1]), default: nx^-2
-  ``lambda_spec`` : Spectral Regularization parameter (weight)
   (Eq.29,31 [1]), default: 1e-2
-  ``lambda_L1`` : L1 constraint, emphasize sparsity of object (as stars
   field), default: 0
-  ``epsilon`` : Ridge/tikhonov parameter epilson \|X\|^2\_2 (Eq.29,31
   [1]), default: 1e-6
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
-  ``FOV`` : Field Of View User parameter, must be in ArcSecond, default: 40 mas
-  ``mask3D``: Support constraint, can be a path to a fits file, an user
   Array, an empty Array is default and returns a unconstraint support,
   or can be initialized with the ``mask(...)`` function. In case of
   array dimensions are checked and corrected, default: no constraint
-  ``xinit3D``: Initial Estimate, as for the mask, this parameter accept
   a path to a fits files. Can also be an Array, or complexe
   visibilities. In case of array dimensions are checked and corrected,
   default: centered dirac att all wavelengths
-  ``indfile``: allow to reduce the set of file present in the
   ``Folder``, default: all files
-  ``indwvl``: allow to restrict the set of wavelength used (wavelength
   index), default: all wavelengths
-  ``CountPlot``: draw plot at each ``CountPlot`` iterations, default:
   10
-  ``admm``: to run admm or not (true or false), false just initialize
   and give structures to run ``painteradmm(...)``, default: true
-  ``paral``: if true perform some calculus in parallel, useful when
   julia use several core (``julia -p nprocs``), default: true

Optimization engine Structures
------------------------------

``OPTOPT`` structure contains:

- ``ls,scl,gat,grt,vt,memsize,mxvl,mxtr,stpmn,stpmx``: related to
  ```OptimPack`` <https://github.com/emmt/OptimPack>`__. default:

.. code:: julia

	ls=OptimPack.MoreThuenteLineSearch(ftol=1e-4,gtol=0.9)
  	scl=OptimPack.SCALING\_OREN\_SPEDICATO
  	gat=1E-6
  	grt=1E-6
  	vt=false
  	memsize=100
  	mxvl=1000
  	mxtr=1000
  	stpmn=1E-20
  	stpmx=1E+20

if parameters are not setted, default value are used. For example,
calling: ``OIDATA,PDATA,OPTOPT =  painter()`` execute the 3D image
reconstruction algorithm from data stored in all \*.oifits files from
folder "OIFITS" located in the ``Painter.jl`` source folder
(``src/OIFITS``). The parameters are setted to default value with no
support contraint, spatial and spectral regularizations, positivity
constraint, the original estimate is a centered dirac at all
wavelengths.

Useful output data
------------------

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
