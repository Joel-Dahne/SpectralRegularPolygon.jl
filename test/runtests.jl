using ArbExtras
using Arblib
using Test

import SpectralRegularPolygons as SRP

@testset "SpectralRegularPolygons" begin
    # Tests for small N
    include("small_N/RegularPolygon.jl")
    include("small_N/enclosure.jl")

    # Tests for large N
    include("large_N/basic_functions.jl")
    include("large_N/basic_integrals.jl")
    include("large_N/fx_div_x.jl")
    include("large_N/TaylorModel.jl")
    include("large_N/polylog.jl")
    include("large_N/section_2.jl")
end
