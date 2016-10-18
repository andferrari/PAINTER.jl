Getting Started
===============

Installation
------------

``PAINTER.jl`` uses the following registered Julia packages:

* `OptimPack.jl <https://github.com/emmt/OptimPack.jl>`_: the Julia interface to `OptimPack <https://github.com/emmt/OptimPack>`_ for solving the phases proximal operator.
* `OIFITS.jl <https://github.com/emmt/OIFITS.jl>`_: Julia support for OI-FITS (optical interferometry data format).
* `NFFT.jl <https://github.com/tknopp/NFFT.jl>`_: Julia implementation of the Non-equidistant Fast Fourier Transform (NFFT).
* `Wavelets.jl <https://github.com/JuliaDSP/Wavelets.jl>`_: A Julia package for fast wavelet transforms.
* `HDF5.jl <https://github.com/timholy/HDF5.jl>`_: for writing JLD ("Julia data") variables.

They will be *automaticaly* installed during ``PAINTER.jl`` installation.
Note that they require development tools included for example for
OSX in ``Command Line Tools`` and for ubuntu in the ``build-essential`` package.

To install ``PAINTER.jl``, type from a Julia session the following commands:

.. code:: julia

  Pkg.update()
  Pkg.add("PAINTER")

and relax!

It is recommended to install ``PyPlot.jl`` to monitor the iterations of the algorithm when the number
of wavelengths is small, e.g. < 30.  See `PyPlot.jl <https://github.com/stevengj/PyPlot.jl>`_ page.

Usage
-----

To load the ``PAINTER.jl`` module, type from a Julia session:

.. code:: julia

    using PAINTER

If ``PyPlot`` is installed, it will be automatically loaded.

Iteration steps of ``PAINTER.jl`` are parallelized.
To use parallel computing, start Julia with ``nprocs`` local process
and load the module on all process:

.. code:: julia

    $ julia -p nprocs
    julia> using PAINTER

Path to ``OptimPack`` options
-----------------------------

The file ``optpckpt.jl`` located in source file of ``PAINTER`` contains all OptimPack parameters for the phases proximal operator.

  * ``ls``, ``scl``, ``gat``, ``grt``, ``vt``, ``memsize``, ``mxvl``, ``mxtr``, ``stpmn``, ``stpmx``. See  `OptimPack <https://github.com/emmt/OptimPack>`_ for details.

  Default values are:

  .. code:: julia

    ls=OptimPack.MoreThuenteLineSearch(ftol=1e-8,gtol=0.95)
    scl=OptimPack.SCALING\_OREN\_SPEDICATO
    gat=0
    grt=0
    vt=false
    memsize=100
    mxvl=1000
    mxtr=1000
    stpmn=1E-20
    stpmx=1E+20
