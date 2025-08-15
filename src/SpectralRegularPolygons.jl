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

import Arblib: Arf, Arb, Acb, ArbSeries, AcbSeries, ArbPoly, AcbPoly, abs_ubound, add_error
import GeometryBasics: Point, Point2
import SpecialFunctions: besselj, besselj0, besselj1, bessely0, bessely1, zeta, gamma

export RegularPolygon, Eigenfunction

include("arb.jl")

include("Polar.jl")
include("RegularPolygon.jl")

include("MPS/Eigenfunction.jl")
include("MPS/sigma.jl")
include("MPS/mps.jl")

include("MPS/enclosing/maximum_boundary.jl")
include("MPS/enclosing/norm.jl")
include("MPS/enclosing/eigenvalue.jl")

include("precomputed_eigenvalues.jl")

include("large_N/special_functions.jl")
include("large_N/basic_functions.jl")
include("large_N/basic_integrals.jl")
include("large_N/tools.jl")
include("large_N/TaylorModel.jl")
include("large_N/polylog.jl")
include("large_N/section_2.jl")
include("large_N/section_2_1_lemmas.jl")
include("large_N/section_2_2_lemmas.jl")
include("large_N/section_2_3_lemmas.jl")

end # module SpectralRegularPolygons
