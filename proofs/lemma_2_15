### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 4cf93e40-85b4-11f0-0b09-df73f8553fa5
begin
    using Pkg
    Pkg.activate("..", io = devnull)
    using SpectralRegularPolygons
    using CairoMakie
    using Arblib
    using ArbExtras
    using OhMyThreads
    using PlutoUI
    using SpecialFunctions

    import SpectralRegularPolygons as SRP
    import ProgressLogging: @withprogress, @logprogress

    setprecision(Arb, 128)
end

# ╔═╡ 33a8bf7d-e4d2-4051-ba4a-16977a14d443
md"""
# Proof of Lemma 2.15
This notebook contains the computer-assisted part of the proof of Lemma 2.15.
"""

# ╔═╡ 5fc3ce43-cfcc-4111-aff2-ce878397f097
md"""
## Goal
We want to prove that for $N \geq N_0$, the regular $N$-sided polygon of area $\pi$ contains a disc of radius $R_{\text{inner}}$ and is contained in a disc of size $R_{\text{outer}}$. Here $N_0$, $R_{\text{inner}}$ and $R_{\text{outer}}$ are given by:
"""

# ╔═╡ a59af2f5-80e6-4b16-9961-d2d62eb6d1be
N₀ = SRP.N₀

# ╔═╡ bb4667b5-bb3c-4cbf-8216-a8f05d8051e3
R_inner = SRP.R_inner

# ╔═╡ 945657ac-d0e5-433a-813a-40bfdf4cbcd4
R_outer = SRP.R_outer

# ╔═╡ 501b94c5-e1c0-43ec-bbc0-74ca2fc6ebf1
md"""
## Proof
The inradius of a regular $N$-sided polygon of area $\pi$ is given by

$$I_N = \sqrt{\frac{\pi}{N\tan(\pi/N)}}$$

and the circumradius by

$$C_N = \sqrt{\frac{\pi}{\frac{N}{2}\sin(2\pi/N)}}.$$

The inradius in an increasing function in $N$ and the circumradius in a decreasing function in $N$. It therefore suffices to verify that

$$I_{N_0} > R_{\text{inner}}$$

and

$$C_{N_0} < R_{\text{outer}},$$

and the statement follows for all $N \geq N_0$. Computing $I_{N_0}$ and $C_{N_0}$ we get:
"""

# ╔═╡ 6ece39a2-2a2f-43b8-b999-0a3ef450a04c
I_N_0 = sqrt(π / (N₀ * tanpi(Arb(1 // N₀))))

# ╔═╡ bdcba7fd-86d5-4a88-ab79-e3f83dcac661
C_N_0 = sqrt(π / (N₀ // 2 * sinpi(Arb(2 // N₀))))

# ╔═╡ bfc91c06-da88-4fdb-afd1-4a8a4d9eeeb1
md"""
Finally, we verify that the required bounds hold:
"""

# ╔═╡ 2bbfb9e5-db35-4248-b3e9-1db3aaaf88c9
I_N_0 >= Arb(R_inner)

# ╔═╡ 789a0a04-5aa5-48de-84a0-277631076b87
C_N_0 <= Arb(R_outer)

# ╔═╡ Cell order:
# ╟─33a8bf7d-e4d2-4051-ba4a-16977a14d443
# ╠═4cf93e40-85b4-11f0-0b09-df73f8553fa5
# ╟─5fc3ce43-cfcc-4111-aff2-ce878397f097
# ╠═a59af2f5-80e6-4b16-9961-d2d62eb6d1be
# ╠═bb4667b5-bb3c-4cbf-8216-a8f05d8051e3
# ╠═945657ac-d0e5-433a-813a-40bfdf4cbcd4
# ╟─501b94c5-e1c0-43ec-bbc0-74ca2fc6ebf1
# ╠═6ece39a2-2a2f-43b8-b999-0a3ef450a04c
# ╠═bdcba7fd-86d5-4a88-ab79-e3f83dcac661
# ╟─bfc91c06-da88-4fdb-afd1-4a8a4d9eeeb1
# ╠═2bbfb9e5-db35-4248-b3e9-1db3aaaf88c9
# ╠═789a0a04-5aa5-48de-84a0-277631076b87
