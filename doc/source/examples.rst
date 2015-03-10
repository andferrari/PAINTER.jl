examples
========

If parameters are not set, default values are used.
For example, calling: ``OIDATA,PDATA,OPTOPT =  painter()`` execute the 3D image
reconstruction algorithm from data stored in all \*.oifits files from
folder "OIFITS" located in the ``Painter.jl`` source folder
(``src/OIFITS``). The parameters are setted to default value with no
support contraint, spatial and spectral regularizations, positivity
constraint, the original estimate is a centered dirac at all
wavelengths.

User parameters and single algorithm execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The folowing parameters are setted by the user, the initial estimate is
the default, drawing is enabled and parallel computing is disabled.
Painter will take oifits informations from all files and will restrict
the analysis on the first 29 wavelength. admm is enabled by default and
will run the algorithm for 100 iterations.

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
    Mynbitermax   = 100
    Mypar         = false

The support constraint is defined by a circle

-  ``Mymask3D = mask(Mynx,int(Mynx/2 -3),"circ")``

The other parameters are setted with default values. Painter is then
executed:

.. code:: julia

    OIDATA,PDATA,OPTOPT = painter(Folder=MyFolder, nbitermax=Mynbitermax, nx=Mynx, lambda_spat=Mylambda_spat=Mylambda_spat, lambda_spec=Mylambda_spec, rho_y= rho_y, rho_spat= rho_spat, rho_spec= rho_spec, rho_ps= rho_ps, alpha= alpha, beta=Mybeta, eps1=Myeps1, eps2=Myeps2, FOV= MyFOV, indfile, indwvl=Myindwvl, paral=Myparal)

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
