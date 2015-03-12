Getting Started
===============

Installation
------------

``Painter.jl`` uses the library
`OptimPack <https://github.com/emmt/OptimPack>`_ for solving the
phases proximal operator. The Julia interface to OptimPack is not yet a
registered Julia package. To install
`OptimPack.jl <https://github.com/emmt/OptimPack.jl>`_ type from a
Julia session the following commands:

.. code:: julia

    Pkg.clone("https://github.com/emmt/OptimPack.jl.git")
    Pkg.build("OptimPack")

``Pkg.build`` requires developpement tools included for example for
OSX in ``Command Line Tools`` and for ubuntu in the ``build-essential`` package.

Installation of PAINTER is then as simple as typing:

.. code:: julia

    Pkg.clone("https://github.com/andferrari/Painter.jl.git")

It is recommended to install ``PyPlot.jl`` to monitor the iterations of the algorithm when the number
of wavelength is small, e.g. < 30.  See `PyPlot.jl <https://github.com/stevengj/PyPlot.jl>`_ page.

Usage
-----

To load the ``Painter.jl`` module, type from a Julia session:

.. code:: julia

    using Painter

If ``PyPlot`` is installed, it will be automatically loaded.

Some iteration steps of ``Painter.jl`` are parallelized.
To use parallel computing, start Julia with ``nprocs`` local process
and load the module on all process:

.. code:: julia

    $ julia -p nprocs
    julia> @everywhere using Painter


Dependencies
------------

``Painter.jl`` uses the following registered Julia package. They will be
automaticaly installed during ``Painter.jl`` installation.

* `FITSIO.jl <https://github.com/JuliaAstro/FITSIO.jl>`_: Julia support for OI-FITS (optical interferometry data format).
* `NFFT.jl <https://github.com/tknopp/NFFT.jl>`_: Julia implementation of the Non-equidistant Fast Fourier Transform (NFFT).
* `Wavelets.jl <https://github.com/JuliaDSP/Wavelets.jl>`_: A Julia package for fast wavelet transforms.
* `HDF5.jl <https://github.com/timholy/HDF5.jl>`_: for writing JLD ("Julia data") variables.
