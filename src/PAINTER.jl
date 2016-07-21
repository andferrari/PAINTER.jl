#
# PAINTER.jl --
#
# Polychromatic opticAl INTErferometric Reconstruction software in Julia.
#
#------------------------------------------------------------------------------
#
# This file is part of PAINTER.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015, Antony Schutz, Andre Ferrari.
#
#------------------------------------------------------------------------------
module PAINTER
    using JLD
    using OIFITS
    using OptimPack
    using NFFT
    using Wavelets

    export painter, mask, paintersave, painterload, painterfitsexport, painterplotfct, WT

    include("paintertype.jl")

    export PAINTER_Input, PAINTER_Data

    include("paintertools.jl")
    include("painterio.jl")
    include("painteroifits.jl")
    include("painterconstmat.jl")
    include("paintercheckinit.jl")
    include("painterplot.jl")
end
