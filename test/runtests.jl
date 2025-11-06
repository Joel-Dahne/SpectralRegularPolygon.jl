using ArbExtras
using Arblib
using Test

import SpectralRegularPolygons as SRP

@testset "SpectralRegularPolygons" begin
    # Tests for small N
    include("small_N/RegularPolygon.jl")

    # Tests for large N
    include("large_N/basic_integrals.jl")
    include("large_N/polylog.jl")
end
