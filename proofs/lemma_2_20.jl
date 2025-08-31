### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ a8ef131c-85bb-11f0-28cf-7b4c670ed13b
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

# ╔═╡ bfb90022-f09a-46a0-b0aa-5d1ffdd7a908
md"""
# Proof of Lemma 2.20
This notebook contains the computer-assisted part of the proof of Lemma 2.20.
"""

# ╔═╡ d2d34a83-67a8-4d71-8b06-171fefb55487
md"""
## Goal
We want to prove that for $N \geq N_0$ we have

$$\frac{\lambda_{\text{app}}(N)}{1 + \varepsilon'(N)} > \frac{\lambda_{\text{app}}(N+1)}{1 - \varepsilon'(N+1)}$$

and

$$\frac{\lambda_{\text{app}}(N)(1-\varepsilon'(N+1))}{(1+\varepsilon'(N))\lambda_{\text{app}}(N+1)} > \frac{\lambda_{\text{app}}(N+1)(1+\varepsilon'(N+2))}{(1-\varepsilon'(N+1))\lambda_{\text{app}}(N+2)},$$

where $N_0$ is given by:
"""

# ╔═╡ 1e814a20-d22d-4c0d-9003-6712a0a5d205
N₀ = SRP.N₀

# ╔═╡ c1de3626-13bd-4a43-ac5c-862c8a8524a4
md"""
## Plots
"""

# ╔═╡ ddaa5aed-4100-4c91-a50b-69a6f7f6aac7
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N^{-1}")
    lines!(ax, range(Arb(0), 1 // N₀, 200), SRP.λ_sup_m_λ_inf)
    fig
end

# ╔═╡ bbd13eae-7074-4618-947f-359a24263951
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N^{-1}")
    lines!(ax, range(Arb(0), 1 // N₀, 200), SRP.q_sup_m_q_inf)
    fig
end

# ╔═╡ 2eadc7be-19a4-4520-b11e-80446cfe0f73
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N^{-1}")
    lines!(
        ax,
        range(Arb(0), 1 // N₀, 200),
        ArbExtras.derivative_function(SRP.λ_sup_m_λ_inf, 4),
    )
    fig
end

# ╔═╡ 8afdf25a-62f4-4dd4-9ec3-f6b409d3422d
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N^{-1}")
    lines!(
        ax,
        range(Arb(0), 1 // N₀, 200),
        ArbExtras.derivative_function(SRP.q_sup_m_q_inf, 5),
    )
    fig
end

# ╔═╡ ecdfb51c-bd88-4d2b-86fd-d6c46340bf28
md"""
## Proof
Let us start by proving that

$$\frac{\lambda_{\text{app}}(N)}{1 + \varepsilon'(N)} > \frac{\lambda_{\text{app}}(N+1)}{1 - \varepsilon'(N+1)}.$$


"""

# ╔═╡ 76af3933-de21-45f0-9310-028eaf31465f
λ_sup_m_λ_inf_expansion = SRP.λ_sup_m_λ_inf(ArbSeries((0, 1), degree = 4))

# ╔═╡ b20f3e75-b475-4cac-b69d-b6611ac82363
iszero(λ_sup_m_λ_inf_expansion[0])

# ╔═╡ 6a2b85e8-510e-499d-be6a-65dffb2934ff
iszero(λ_sup_m_λ_inf_expansion[1])

# ╔═╡ 5d9a1d3b-b2b6-4e5b-a769-5d35e3e64963
iszero(λ_sup_m_λ_inf_expansion[2])

# ╔═╡ 0d216fba-d53a-47ba-886f-dc8902b0b0c2
iszero(λ_sup_m_λ_inf_expansion[3]) # PROVE

# ╔═╡ f7c11696-e359-48c6-bf29-b799b3a26957
# ╠═╡ disabled = true
#=╠═╡
λ_sup_m_λ_inf_d4 = ArbExtras.minimum_enclosure(
	ArbExtras.derivative_function(SRP.λ_sup_m_λ_inf, 4),
	Arf(0),
	ubound(Arb(1 // N₀)),
	verbose = true,
)
  ╠═╡ =#

# ╔═╡ 1b5615d5-ff94-41c6-8123-d11bdae8af08
#=╠═╡
Arblib.ispositive(λ_sup_m_λ_inf_d4)
  ╠═╡ =#

# ╔═╡ 4371e67c-6c2e-461c-bbb5-afcb8ad7f142
md"""
Next we prove that

$$\frac{\lambda_{\text{app}}(N)(1-\varepsilon'(N+1))}{(1+\varepsilon'(N))\lambda_{\text{app}}(N+1)} > \frac{\lambda_{\text{app}}(N+1)(1+\varepsilon'(N+2))}{(1-\varepsilon'(N+1))\lambda_{\text{app}}(N+2)}.$$
"""

# ╔═╡ c427a061-0bf9-4e40-947e-2bd093ad686a
q_sup_m_q_inf_expansion = SRP.q_sup_m_q_inf(ArbSeries((0, 1), degree = 5))

# ╔═╡ 40f8bc3e-9d69-4e95-800c-7f33e5187b8c
iszero(q_sup_m_q_inf_expansion[0])

# ╔═╡ 218d06e7-c755-41b3-bfaa-1dfd6a3ccded
iszero(q_sup_m_q_inf_expansion[1])

# ╔═╡ f167c7ae-09bb-42c7-952b-9b6a24b165ee
iszero(q_sup_m_q_inf_expansion[2])

# ╔═╡ 15efde20-171e-4283-9425-6012e3e1036b
iszero(q_sup_m_q_inf_expansion[3]) # PROVE

# ╔═╡ 64c4498a-979b-4cad-bfb9-b193ca2c93d1
iszero(q_sup_m_q_inf_expansion[4]) # PROVE

# ╔═╡ dbf0f411-344f-4d6c-b031-b5c61022cbd4
a = Arf(1e-3)

# ╔═╡ a14c2a9a-4561-4420-b7aa-fe386bc84ea2
q_sup_m_q_inf_d4_0_a = ArbExtras.minimum_enclosure(
    ArbExtras.derivative_function(SRP.q_sup_m_q_inf, 5),
    Arf(0),
    a,
    verbose = true,
)

# ╔═╡ 5169dc6c-8189-460a-a7c4-37bd4360f7d3
Arblib.ispositive(q_sup_m_q_inf_d4_0_a)

# ╔═╡ 94b731fb-019c-43e1-b0b3-519d0f801529
q_sup_m_q_inf_a_inv_N =
    ArbExtras.minimum_enclosure(SRP.q_sup_m_q_inf, a, ubound(Arb(1 // N₀)), verbose = true)

# ╔═╡ cdeb73d6-600e-49b4-8dd9-2329d4b9dc9d
Arblib.ispositive(q_sup_m_q_inf_a_inv_N)

# ╔═╡ Cell order:
# ╟─bfb90022-f09a-46a0-b0aa-5d1ffdd7a908
# ╠═a8ef131c-85bb-11f0-28cf-7b4c670ed13b
# ╟─d2d34a83-67a8-4d71-8b06-171fefb55487
# ╠═1e814a20-d22d-4c0d-9003-6712a0a5d205
# ╟─c1de3626-13bd-4a43-ac5c-862c8a8524a4
# ╟─ddaa5aed-4100-4c91-a50b-69a6f7f6aac7
# ╟─bbd13eae-7074-4618-947f-359a24263951
# ╟─2eadc7be-19a4-4520-b11e-80446cfe0f73
# ╟─8afdf25a-62f4-4dd4-9ec3-f6b409d3422d
# ╠═ecdfb51c-bd88-4d2b-86fd-d6c46340bf28
# ╠═76af3933-de21-45f0-9310-028eaf31465f
# ╠═b20f3e75-b475-4cac-b69d-b6611ac82363
# ╠═6a2b85e8-510e-499d-be6a-65dffb2934ff
# ╠═5d9a1d3b-b2b6-4e5b-a769-5d35e3e64963
# ╠═0d216fba-d53a-47ba-886f-dc8902b0b0c2
# ╠═f7c11696-e359-48c6-bf29-b799b3a26957
# ╠═1b5615d5-ff94-41c6-8123-d11bdae8af08
# ╟─4371e67c-6c2e-461c-bbb5-afcb8ad7f142
# ╠═c427a061-0bf9-4e40-947e-2bd093ad686a
# ╠═40f8bc3e-9d69-4e95-800c-7f33e5187b8c
# ╠═218d06e7-c755-41b3-bfaa-1dfd6a3ccded
# ╠═f167c7ae-09bb-42c7-952b-9b6a24b165ee
# ╠═15efde20-171e-4283-9425-6012e3e1036b
# ╠═64c4498a-979b-4cad-bfb9-b193ca2c93d1
# ╠═dbf0f411-344f-4d6c-b031-b5c61022cbd4
# ╠═a14c2a9a-4561-4420-b7aa-fe386bc84ea2
# ╠═5169dc6c-8189-460a-a7c4-37bd4360f7d3
# ╠═94b731fb-019c-43e1-b0b3-519d0f801529
# ╠═cdeb73d6-600e-49b4-8dd9-2329d4b9dc9d
