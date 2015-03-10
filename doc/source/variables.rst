Variables and structures
========================

Variables
---------

  -  ``nbitermax``: number of ADMM iterations. Default ``1000``
  -  ``aff``: if ``aff=true`` plots are enabled using ``PyPlot.jl``. Default: ``false``.

Variables in ``OIDATA`` structure
----------------------------------

``OIDATA``: contains all oifits information and user defined parameters.

**Execution Variables:**

  - ``paral``: if ``paral=true`` the ADMM step which reconstructs the object for each wavelength is computed in parallel.
    In this case ``julia``must be started with

    .. code:: bash

      $ julia -p nprocs

    where ``nprocs``denotes the number of processes. Default: true
  - ``admm``: if ``admm=false`` the function only initializes the structures. The function ``painter`` can be used after to iterate
    the ADMM algorithm. Default: ``true``.
  -  ``CountPlot``: draw plot at each ``CountPlot`` iterations. Default: 10

**Data related variables:**

  - ``Folder``: path to the folder containing oifits/fits files. Default: ``src/OIFITS``.
  - ``indfile``: allows to chose the set of oifits/fits files processed in ``Folder``. ``indfile`` is an Array of Int containnig the alphabetical order of the files. Default: all files.
  -  ``indwvl``: allowd to restrict the set of processed wavelengths. ``indwvl`` is an  Array of Int containnig the order of the wavelengths. Default: all wavelengths.
  - ``nx``: image size in pixels (image of size nx\ :sup:`2`). Default: ``64``.
  - ``FOV``: Field Of View of the reconstructed image in ArcSecond. Default: ``40e-3``
  - ``mask3D``: Binary mask defining the support constraint. ``mask3D`` can be

    - a path to a fits file,
    - an Array,
    - an empty Array (no constraint).

    ``mask3D`` can be set by function ``mask``. Default: no constraint.

  -  ``xinit3D``: Initial Estimate of the object or of the complex visibilities. ``xinit3D`` can be:

    - a path to a fits files containing the object,
    - an Array containing the object,
    - and Array containing the complex visibilities

    Default: centered dirac for all wavelengths.


**ADMM algorithm parameters:**

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


Variables in ``OPTOPT`` structures
----------------------------------

``OPTOPT``: contains all OptimPack parameters for the phases proximal operator.


  - ``ls,scl,gat,grt,vt,memsize,mxvl,mxtr,stpmn,stpmx``: related to
  `OptimPack <https://github.com/emmt/OptimPack>`_.

  Default:

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


Useful output variables in ``PDATA`` structure
----------------------------------------------

Useful Array in ``PDATA`` are

  - ``PDATA.x``: reconstruced 3D images
  - ``PDATA.w``: positivity and support contraint. These constraints can be applied to ``PDATA.x``
    with ``PDATA.x.*(PDATA.w.>0)``.
  - ``PDATA.Fx``: non uniform Fourier transform of the reconstructed 3D images.
  - ``PDATA.crit1``: the primal residual of the ADMM algorithm.
  - ``PDATA.crit2``: the dual residual of the ADMM algorithm.
  - ``PDATA.ind``: number of iterations, useful to re-run algorithm.
