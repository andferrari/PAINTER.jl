Variables and structures
========================

If parameters are not set, default values are used. For example

.. code:: julia

  OIDATA,PDATA =  painter()

calls ``painter`` with all default values.

Variables
---------

These two variables cannot be included in ``OIDATA`` structure.

* ``nbitermax``: number of ADMM iterations. Default: ``1000``.
* ``aff``: if ``aff=true`` plots are enabled using ``PyPlot.jl``. Default: ``false``.

Variables in ``OIDATA`` structure
----------------------------------

The structure ``OIDATA`` contains all OIFITS information and user defined parameters.

**Execution Variables:**

* ``admm``: if ``admm=false`` the function only initializes the structures. The function ``painter`` can be used after to iterate the ADMM algorithm. Default: ``true``.
* ``CountPlot``: draw plot at each ``CountPlot`` iterations. Default: ``10``.
* ``PlotFct``: is a user defined function which is called at each ``CountPlot`` iterations. This function must respect the input argument of ``painterplotfct`` function and must call ``PyPlot``, see :ref:`examples-label`  section. Default: ``painterplotfct``.

**Initialization and initial estimate:**

* ``autoinit``: if ``autoinit=true`` some parameters are automatically set or rescaled. Default: ``true``.

The parameters which are automaticaly initialized are ``alpha``, ``beta``, ``rho_y_xi`` and ``rho_y_gamma``.
They corresponds to parameters related to proximal operator for squared visibilities and for phases differences.
Regularization parameters ``lambda_spat`` and ``lambda_spec`` are rescaled to be invariant with user parameters.
``lambda_spat`` is divided by the number of pixels and ``lambda_spec`` by the number of wavelength.
The total flux is also normalized to allow the use of almost predefined parameters.
The initial estimate is rescaled by the flux of the data.

**Data and image related variables:**

* ``Folder``: path to the folder containing OIFITS/FITS files. Default: ``./OIFITS``. If ``./OIFITS`` does not exists ``src/OIFITS`` in ``PAINTER.jl/`` default installation folder, containing FITS files for the demo, is used.
* ``indfile``: allows to chose the set of OIFITS/FITS files in ``Folder`` that will be processed. ``indfile`` is an ``Array{Int64,1}`` containing the indexes of the files in alphabetical order. Default: all files.
* ``indwvl``: allows to chose the set of processed wavelengths. ``indwvl`` is an ``Array{Int64,1}`` containing the indexes of the wavelengths in increasing order. Default: all wavelengths.
* ``nx``: image size in pixels (the size of the image is nx\ :sup:`2`). Default: ``64``.
* ``FOV``: Field Of View of the reconstructed image in ArcSecond. Default: ``40e-3``.
* ``mask3D``: Binary mask defining the support constraint. ``mask3D`` can be:

  * a path to a FITS file,
  * an Array,
  * an empty Array (no constraint).

  ``mask3D`` can be set by the function ``mask``. Default: no constraint.

* ``xinit3D``: Initial estimate of the object or of the complex visibilities. ``xinit3D`` can be:

  * a path to a FITS files containing the object,
  * an Array containing the object,
  * and Array containing the complex visibilities.

  Default: centered Dirac functions at all wavelengths.

* ``dptype`` define the kind of matrix difference used to generate differential phase, can be parameterized by ``dpprm`` :

  * ``"all"`` the difference between the first wavelength and all others (1-2, 1-3, ...), see  Eqs. 35
  * ``"diag"`` the difference between all consecutive wavelengths (1-2, 2-3, ...)
  * ``"ref"`` the same as ``"all"`` but with a reference channel defined by ``dpprm``, the same as ``"all"`` if ``dpprm``=1
  * ``"frame"`` the difference between wavelength are performed inside non overlapping window with a size ``dpprm``
  * ``"sliding"`` the difference between wavelength are performed using a sliding window with a size ``dpprm``

  Default: if not given the default matrix difference is ``"all"``, for details about other methods see [3].

**ADMM algorithm parameters:**

* ``alpha``: weight for squared visibilities modulus data fidelity term, see Eqs. 25, 31 in [1]_. Default: ``1``.
* ``beta``: weight for phases (closures and differential) data fidelity term, see Eqs. 25,31 in [1]_. Default: ``1``.
* ``lambda_spat``: Spatial regularization parameter, see Eqs. 29, 31 in [1]_. Default: nx\ :sup:`-2`.
* ``lambda_spec``: Spectral regularization parameter, see Eqs. 29, 31 in [1]_. Default: ``1e-2``.

* ``rho_y``: ADMM parameter for data fidelity,see  Eqs. 35, 50-52 in [1]_. Default: ``10``.
* ``rho_spat``: ADMM parameter for Spatial regularization, see Eqs. 25, 31 in [1]_. Default: ``1``, (``0`` to disable).
* ``rho_spec``: ADMM parameter for Spectral regularization, see Eqs. 42, 55 in [1]_. Default: ``1``, (``0`` to disable).
* ``rho_ps``: ADMM parameter for positivity constraint, see Eq. 47, 54 in [1]_. Default: ``1``, (``0`` to disable).

Secondary or specific paramaters:
  The defaults values of these parameteres are tuned for the general cases. Nevertheless, the user may modified them for specific applications.

  * ``lambda_L1``: regularization parameter for an l\ :sub:`1` constraint on the image. l\ :sub:`1` constraint emphasizes sparsity of objects (e.g. stars field). Default: ``0``.
  * ``Wvlt``: array of wavelets basis for spatial regularization, see [2]_.  See `Wavelets.jl <https://github.com/JuliaDSP/Wavelets.jl>`_ for definitions. Default: first 8 Daubechies wavelets and Haar wavelets. ``Wvlt = [WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar]``.
  * ``epsilon``: Ridge/Tikhonov regularization parameter, see Eqs. 29, 31 in [1]_. Default: ``1e-6``.
  * ``eps1``: stopping criterium  for primal residual  in ADMM algorithm. Default: ``1e-6``.
  * ``eps2``: stopping criterium for dual residual in ADMM algorithm. Default: ``1e-6``.


Constant in ``OIDATA`` structure
--------------------------------

The structure ``OIDATA``: contains also constants related to the data and
extracted from OIFITS files.

* ``nb``: number of bases.
* ``nw``: number of wavelength.
* ``U``: the U spatial frequencies matrix.
* ``V``: the V spatial frequencies matrix.
* ``P``: squared visibilities Matrix.
* ``W``: squared visibilities variance Matrix.
* ``T3``: phases closure matrix.
* ``T3err``: phases closure variance matrix.
* ``DP``: differential phases vector.
* ``DPerr``: differential phases variance vector.
* ``Xi``: dictionary of phases difference Vector.
* ``K``: dictionary of phases difference variance vector.

For matrices, the column index is associated to the wavelength index.

Variables in ``PDATA`` structure
--------------------------------

Useful outputs in the structure ``PDATA`` are:

* ``PDATA.x``: the reconstructed 3D images !
* ``PDATA.w``: positivity and support constraint. These constraints can be applied to ``PDATA.x`` with ``PDATA.x.*(PDATA.w.>0)``.
* ``PDATA.Fx``: non uniform Fourier transform of the reconstructed 3D images.
* ``PDATA.H``: dictionary of phases to phases differences sparse matrix.
* ``PDATA.crit1``: the primal residual of the ADMM algorithm.
* ``PDATA.crit2``: the dual residual of the ADMM algorithm.
* ``PDATA.ind``: number of iterations, useful to re-run algorithm.

References
----------

.. [1] Schutz, A., Ferrari, A., Mary, D. Soulez, F., Thiébaut, E., Vannier, M. "PAINTER: a spatio-spectral image reconstruction algorithm for optical interferometry". JOSA A. Vol. 31, Iss. 11, pp. 2356–2361, (2014). `arXiv <http://arxiv.org/abs/1407.1885>`_
.. [2] Schutz, A., Ferrari, A., Mary, D., Thiébaut, E., Soulez, F. "Large scale 3D image reconstruction in optical interferometry". EUSIPCO, 2015, Nice. `arXiv <http://arxiv.org/abs/1503.01565>`_
.. [3] Schutz, A., Ferrari, A., Thiébaut, E., Soulez, F., Vannier, M., Mary D. "Interbands phase models for polychromatic image reconstruction in optical interferometry". SPIE, 2016, Edinburgh.
