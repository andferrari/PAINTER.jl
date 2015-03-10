Variables and structures
========================

Two parameters are independent of structures:

-  ``nbitermax``: number of
   ADMM
   iteration, default 1000
-  ``aff``: true for plot drawing using ``PyPlot.jl``, default: false

``OIDATA`` structure
--------------------

  ``OIDATA`` structure contains the following fields:

  - ``Folder``: path to the folder containing oifits/fits files. Default: ``src/OIFITS``.
  - ``nx``: image size in pixels (image of size nx\ :sup:`2`). Default: ``64``.
  - ``Wvlt``: list of wavelets basis. See `Wavelets.jl <https://github.com/JuliaDSP/Wavelets.jl>`_. Default: first 8 Daubechies wavelets and Haar wavelets.
  - ``lambda_spat``: Spatial regularization parameter (weight) (Eqs. 29, 31 in ref. [1]). Default: nx\ :sup:`-2`.
  - ``lambda_spec``: Spectral regularization parameter (weight) (Eqs. 29, 31 in ref. [1]). Default: ``1e-2``.
  - ``lambda_L1``: l1 regularization parameter (weight). l1 constraint emphasizes sparsity of objects (e.g. stars field). Default: ``0``.
  - ``epsilon``: Ridge/tikhonov regularization parameter (Eqs. 29, 31 in ref. [1]). Default: ``1e-6``.
  - ``rho_y``: ADMM parameter for data fidelity (convergence rate) (Eq.35, 50-52 in ref. [1]). Default: ``1``.
  - ``rho_spat``: ADMM parameter for Spatial regularization (Eqs. 25, 31 in ref. [1]). Default: ``1``.
  - ``rho_spec``: ADMM parameter for Spectral regularization (Eqs. 42, 55 in ref. [1]). Default: ``1``.
  - ``rho_ps``: ADMM parameter for positivity constraint (Eq. 47, 54 in ref. [1]). Default: ``1``.
  - ``alpha``: weight of absolute squared visibilities data fidelity term (Eqs. 25, 31 in ref. [1]). Default: ``1``.
  - ``beta``: weight for phases (closures and differential) data fidelity term (Eqs. 25,31 in ref. [1]). Default: ``1``.
  - ``eps1``: Primal Residual stopping criterium in ADMM algorithm. Default: ``1e-6``.
  - ``eps2``: Dual Residual stopping criterium in ADMM algorithm. Default: ``1e-6``.
  - ``FOV``: Field Of View of the reconstructed image in ArcSecond. Default: ``40e-3``
  - ``mask3D``: Binary mask defining the support constraint. ``mask3D`` can be

    - a path to a fits file,
    - an user Array,
    - an empty Array (no constraint).

    ``mask3D`` can be set by function ``mask(...)``. Default: no constraint.
  -  ``xinit3D``: Initial Estimate of the object or of the complex visibilities. ``xinit3D`` can be

  , as for the mask, this parameter accept a path to a fits files. Can also be an Array, or complexe
    visibilities. In case of array dimensions are checked and corrected, default: centered dirac att all wavelengths

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
   ADMM
   algorithm
-  ``PDATA.ind`` the iteration indice, useful to re-run algorithm
