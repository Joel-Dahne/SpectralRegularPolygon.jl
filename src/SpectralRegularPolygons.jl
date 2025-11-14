module SpectralRegularPolygons

import Arblib
import ArbExtras
import CSV
import DataFrames
import GenericLinearAlgebra
import GeometryBasics
import LaTeXStrings
import LinearAlgebra
import OhMyThreads
import Optim
import Random

import Arblib:
    Arf,
    Arb,
    Acb,
    ArbSeries,
    AcbSeries,
    ArbPoly,
    AcbPoly,
    radius,
    midpoint,
    ubound,
    abs_ubound,
    add_error,
    getinterval
import GeometryBasics: Point, Point2
import SpecialFunctions: besselj, besselj0, besselj1, bessely0, bessely1, zeta, gamma

# This is the point at which we switch between the case for small and
# large N.
const N₀ = 64

include("assert_proof.jl")

# Some generic code use for both small and large N
include("arb.jl")

# Code for small N
include("small_N/Polar.jl")
include("small_N/RegularPolygon.jl")
include("small_N/MPS/Eigenfunction.jl")
include("small_N/MPS/sigma.jl")
include("small_N/MPS/mps.jl")
include("small_N/MPS/enclosing/maximum_boundary.jl")
include("small_N/MPS/enclosing/norm.jl")
include("small_N/MPS/enclosing/eigenvalue.jl")
include("small_N/precomputed_eigenvalues.jl")

# Code for large N
include("large_N/constants.jl")
include("large_N/special_functions.jl")
include("large_N/basic_functions.jl")
include("large_N/basic_integrals.jl")
include("large_N/fx_div_x.jl")
include("large_N/TaylorModel.jl")
include("large_N/polylog.jl")
include("large_N/section_2.jl")
include("large_N/V.jl")
include("large_N/lemma_2_10.jl")

end # module SpectralRegularPolygons
