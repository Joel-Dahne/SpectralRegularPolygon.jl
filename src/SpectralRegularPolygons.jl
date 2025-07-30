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

import Arblib: Arf, Arb, Acb, ArbSeries, AcbSeries, AcbPoly
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

include("large_N/tools.jl")
include("large_N/TaylorModel.jl")
include("large_N/polylog.jl")
include("large_N/section_2.1.jl")
include("large_N/section_2.1_lemmas.jl")

end # module SpectralRegularPolygons
