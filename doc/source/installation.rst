Getting Started
===============

Installation
------------

``Painter.jl`` uses the library
`OptimPack <https://github.com/emmt/OptimPack>`_ for solving the
phases proximal operator. The Julia interface to OptimPack is not a
registered Julia package. To install
`OptimPack.jl <https://github.com/emmt/OptimPack.jl>`_ type from a
Julia session the following commands:

.. code:: julia

    julia> Pkg.clone("https://github.com/emmt/OptimPack.jl.git")
    julia> Pkg.build("OptimPack")

Installation of PAINTER is then as simple as typing:

.. code:: julia

    julia> Pkg.clone("https://github.com/andferrari/Painter.jl.git")

It is recommended to install ``PyPlot.jl`` to monitor the iterations of the algorithm when the number
of wavelength is small, e.g. < 30.  See `PyPlot.jl <https://github.com/stevengj/PyPlot.jl>`_ page.

Dependencies
------------

``Painter.jl`` uses the following registered Julia package. They will be
automaticaly installed during ``Painter.jl`` installation.

- `FITSIO.jl <https://github.com/JuliaAstro/FITSIO.jl>`_: Julia
   support for OI-FITS (optical interferometry data format).
- `NFFT.jl <https://github.com/tknopp/NFFT.jl>`_: Julia
   implementation of the Non-equidistant Fast Fourier Transform (NFFT).
- `Wavelets.jl <https://github.com/JuliaDSP/Wavelets.jl>`_: A Julia
   package for fast wavelet transforms.
- `HDF5.jl <https://github.com/timholy/HDF5.jl>`_: for writing JLD
   ("Julia data") variables.
