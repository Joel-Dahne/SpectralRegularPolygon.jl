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

import Arblib: Arf, Arb, ArbSeries
import GeometryBasics: Point, Point2
import SpecialFunctions: besselj, bessely0, bessely1

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

end # module SpectralRegularPolygons
