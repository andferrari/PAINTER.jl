Credits
=======

``Painter.jl`` is a julia implementation of PAINTER: Polychromatic
opticAl INTErferometric Reconstruction software described in [1] and
[2].

PAINTER was developped at Laboratoire J.-L. Lagrange, Université de Nice
Sophia, CNRS, Observatoire de la Côte d'Azur, by `Antony
Schutz <http://www.antonyschutz.com>`__ and `André
Ferrari <https://www-n.oca.eu/aferrari>`__.

References
----------

The PAINTER algorithm is described in [1]. The original MATLAB code is
available `here <https://www-n.oca.eu/aferrari/painter/>`__ but the use
of ```Painter.jl`` <https://github.com/andferrari/Painter.jl>`__ is
highly recommended.
```Painter.jl`` <https://github.com/andferrari/Painter.jl>`__ implements
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

.. |Build Status| image:: https://travis-ci.org/andferrari/Painter.jl.svg?branch=master
   :target: https://travis-ci.org/andferrari/Painter.jl
