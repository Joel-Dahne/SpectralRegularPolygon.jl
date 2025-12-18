### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 0ef85c98-edfb-4f07-87b3-77eca945ee78
begin
    using Pkg
    Pkg.activate("..", io = devnull)
    using SpectralRegularPolygons
    using Arblib
    using ArbExtras
    using SpecialFunctions

    import SpectralRegularPolygons as SRP

    setprecision(Arb, 128)
end

# ╔═╡ a5dca01c-dbb8-11f0-3370-f7e17998c69e
md"""
# Proof of Corollary 2.18
"""

# ╔═╡ b00fef7d-6a12-424e-94c4-434df25a7895
md"""
## Goal
We want to prove that for $N \geq N_0$ we have the inequality

$$\sqrt{\lambda_{app}} R_{\text{inner}} < j_{0,1}$$,

where $N_0$ and $R_{\text{unner}}$ are given by
"""

# ╔═╡ 27eb68a5-556f-4a0e-a41d-4926f8e6fb9e
N₀ = SRP.N₀

# ╔═╡ 9b1382f1-f0ad-4439-a64e-6758df7ed3e9
R_inner = SRP.R_inner

# ╔═╡ d649d812-a8be-48c7-8a98-048ab6ba4847
md"""
## Proof
This is a straightforward computation. For the left hand side we get the enclosure
"""

# ╔═╡ 536945f2-2743-4a39-b0c2-c1428dcbde8f
lhs = sqrt(SRP.λ_app(Arb((0, 1 // N₀)))) * Arb(R_inner)

# ╔═╡ d462e6ef-f546-40ca-b709-3c485523d512
md"""
For the right hand side we use that $j_{0,1}$ is the square root of the first eigenvalue of the unit disc, giving us the enclosure
"""

# ╔═╡ 9dbde9ee-1770-4df8-9b52-1d1bf0293371
rhs = sqrt(SRP.λ_disc())

# ╔═╡ a3518d56-2451-4604-aea6-5ea3a96b23c9
md"""
We can now directly verify the condition.
"""

# ╔═╡ 74e20245-b71c-43c6-bd5e-f4f6dfd30755
@assert_proof lhs < rhs

# ╔═╡ ba4f4b4b-6f0b-419f-92a2-2ae9dfb26d15
string(rhs, digits = 10) # For inclusing in the paper we print fewer digits

# ╔═╡ Cell order:
# ╟─a5dca01c-dbb8-11f0-3370-f7e17998c69e
# ╠═0ef85c98-edfb-4f07-87b3-77eca945ee78
# ╟─b00fef7d-6a12-424e-94c4-434df25a7895
# ╠═27eb68a5-556f-4a0e-a41d-4926f8e6fb9e
# ╠═9b1382f1-f0ad-4439-a64e-6758df7ed3e9
# ╟─d649d812-a8be-48c7-8a98-048ab6ba4847
# ╠═536945f2-2743-4a39-b0c2-c1428dcbde8f
# ╟─d462e6ef-f546-40ca-b709-3c485523d512
# ╠═9dbde9ee-1770-4df8-9b52-1d1bf0293371
# ╟─a3518d56-2451-4604-aea6-5ea3a96b23c9
# ╠═74e20245-b71c-43c6-bd5e-f4f6dfd30755
# ╠═ba4f4b4b-6f0b-419f-92a2-2ae9dfb26d15
