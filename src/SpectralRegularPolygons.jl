module SpectralRegularPolygons

import GLMakie
import GenericLinearAlgebra
import GeometryBasics
import LaTeXStrings
import LinearAlgebra
import Optim
import Random

import GeometryBasics: Point, Point2
import SpecialFunctions: besselj

include("Polar.jl")
include("RegularPolygon.jl")

include("MPS/Eigenfunction.jl")
include("MPS/sigma.jl")
include("MPS/mps.jl")
include("MPS/enclosing/maximum_boundary.jl")
include("MPS/enclosing/norm.jl")

include("plotting.jl")

end # module SpectralRegularPolygons
