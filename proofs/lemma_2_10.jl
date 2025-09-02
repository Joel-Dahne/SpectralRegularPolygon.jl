### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ fca89648-853b-11f0-3116-5b610b900eb4
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

# ╔═╡ 18d792d6-8c0c-4214-bd98-670611da3e1e
md"""
# Proof of Lemma 2.10
This notebook contains the computer-assisted part of the proof of Lemma 2.10.
"""

# ╔═╡ d22dc5aa-bac2-4f8a-968b-ece39e295129
N₀ = SRP.N₀

# ╔═╡ bdc6dcff-50a0-49d7-98a2-00b24f5ecfe6
md"""
We want to prove that for $|z| = 1$ we have the following bounds:

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} d_1(z, t)V_4(t)\ dt\right| \leq C_{I,1,4},$$

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} d_2(z, t)V_3(t)\ dt\right| \leq C_{I,2,3},$$

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} d_3(z, t)V_2(t)\ dt\right| \leq C_{I,3,2}$$

and for $N \geq$ $N₀ have

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} K_4(z, t)V_1(t)\ dt\right| \leq C_{I,4,1}.$$

Here $N_0$ and $C_{I,k,l}$ are given by
"""

# ╔═╡ 5c5f73da-4135-412f-a31f-da02450f235a
C_I_1_4 = SRP.C_I_1_4

# ╔═╡ 4b108a4d-cd35-4f5c-ad2f-523d6c3c73d9
C_I_2_3 = SRP.C_I_2_3

# ╔═╡ 11392f25-b796-4f5d-a970-c0a2f82bdcd6
C_I_3_2 = SRP.C_I_3_2

# ╔═╡ 278f9795-6b9b-4780-a0c7-87a12efa8e09
C_I_4_1 = SRP.C_I_4_1

# ╔═╡ b6096ef5-2ce3-4519-97db-6ccb566037f0
md"""
## Plots
"""

# ╔═╡ ef2ecb67-904b-4488-8e47-b5e2fc46eaab
θs_div_πs = range(Arb(0), 2, 100)

# ╔═╡ d24e7344-a2cd-4ba3-a5cf-1290fe285cfc
θs = π * θs_div_πs

# ╔═╡ cf4e04c7-26f2-400b-a885-c3e495e42f78
I_d_1_V_4 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_d_k_V_l(1, 4, SRP.exppii(θ_div_π))
end

# ╔═╡ 67e1ad91-7da5-4985-b8c3-ca1cc84a47fb
I_d_2_V_3 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_d_k_V_l(2, 3, SRP.exppii(θ_div_π))
end

# ╔═╡ a820d24b-9387-4ef5-a959-a4baaf93fdb2
I_d_3_V_2 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_d_k_V_l(3, 2, SRP.exppii(θ_div_π))
end

# ╔═╡ a0bcac27-bbf7-481e-96bd-f1b619e7a54a
I_K_4_V_1 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_K_4_V_1(N₀, SRP.exppii(θ_div_π))
end

# ╔═╡ 7bc1266a-4f01-439b-81fb-5b1e86612ef2
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta")
    band!(ax, θs, lbound.(real(I_d_1_V_4)), ubound.(real(I_d_1_V_4)))
    hlines!(ax, [-Arb(C_I_1_4), Arb(C_I_1_4)])
    fig
end

# ╔═╡ 6187163b-832d-4c9f-862f-a3483505fda6
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta")
    band!(ax, θs, lbound.(real(I_d_2_V_3)), ubound.(real(I_d_2_V_3)))
    hlines!(ax, [-Arb(C_I_2_3), Arb(C_I_2_3)])
    fig
end

# ╔═╡ a7509228-24ea-4395-8bdc-57457c63e28b
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta")
    band!(ax, θs, lbound.(real(I_d_3_V_2)), ubound.(real(I_d_3_V_2)))
    hlines!(ax, [-Arb(C_I_3_2), Arb(C_I_3_2)])
    fig
end

# ╔═╡ bab7bd70-2f53-4421-be57-c79c469c6a5a
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta")
    band!(ax, θs, lbound.(real(I_K_4_V_1)), ubound.(real(I_K_4_V_1)))
    hlines!(ax, [-Arb(C_I_4_1), Arb(C_I_4_1)])
    fig
end

# ╔═╡ 93cf97f7-54c1-4101-93b8-09890db8de86
md"""
## Proof
"""

# ╔═╡ f293357f-b59f-46e8-9b68-a83f1a7cbe42
@time ArbExtras.maximum_enclosure(
	Arf(0),
	Arf(1),
	degree = -1,
	rtol = 1e-3,
	ubound_tol = Arb(C_I_1_4),
	depth_start = 4,
	depth = 30,
	abs_value = true,
	threaded = true,
	verbose = true,
) do θ_div_π
	real(SRP.integral_d_k_V_l(1, 4, SRP.exppii(θ_div_π)))
end

# ╔═╡ 0cfdcd21-461a-4626-b7c8-f6257b9ff3e6
@time ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    ubound_tol = Arb(C_I_2_3),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_d_k_V_l(2, 3, SRP.exppii(θ_div_π)))
end

# ╔═╡ 3d283f17-149d-4e1d-a6d3-0ed64a449490
@time ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    ubound_tol = Arb(C_I_3_2),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_d_k_V_l(3, 2, SRP.exppii(θ_div_π)))
end

# ╔═╡ 67c63d90-496a-4fca-9146-b2ba24ebb4be
@time ArbExtras.maximum_enclosure(
	Arf(0),
	Arf(1),
	degree = -1,
	ubound_tol = Arb(C_I_4_1),
	depth_start = 4,
	depth = 30,
	abs_value = true,
	threaded = true,
	verbose = true,
) do θ_div_π
	real(SRP.integral_K_4_V_1(N₀, SRP.exppii(θ_div_π)))
end

# ╔═╡ Cell order:
# ╟─18d792d6-8c0c-4214-bd98-670611da3e1e
# ╠═fca89648-853b-11f0-3116-5b610b900eb4
# ╟─bdc6dcff-50a0-49d7-98a2-00b24f5ecfe6
# ╠═d22dc5aa-bac2-4f8a-968b-ece39e295129
# ╠═5c5f73da-4135-412f-a31f-da02450f235a
# ╠═4b108a4d-cd35-4f5c-ad2f-523d6c3c73d9
# ╠═11392f25-b796-4f5d-a970-c0a2f82bdcd6
# ╠═278f9795-6b9b-4780-a0c7-87a12efa8e09
# ╟─b6096ef5-2ce3-4519-97db-6ccb566037f0
# ╠═ef2ecb67-904b-4488-8e47-b5e2fc46eaab
# ╠═d24e7344-a2cd-4ba3-a5cf-1290fe285cfc
# ╠═cf4e04c7-26f2-400b-a885-c3e495e42f78
# ╠═67e1ad91-7da5-4985-b8c3-ca1cc84a47fb
# ╠═a820d24b-9387-4ef5-a959-a4baaf93fdb2
# ╠═a0bcac27-bbf7-481e-96bd-f1b619e7a54a
# ╟─7bc1266a-4f01-439b-81fb-5b1e86612ef2
# ╟─6187163b-832d-4c9f-862f-a3483505fda6
# ╟─a7509228-24ea-4395-8bdc-57457c63e28b
# ╟─bab7bd70-2f53-4421-be57-c79c469c6a5a
# ╟─93cf97f7-54c1-4101-93b8-09890db8de86
# ╠═f293357f-b59f-46e8-9b68-a83f1a7cbe42
# ╠═0cfdcd21-461a-4626-b7c8-f6257b9ff3e6
# ╠═3d283f17-149d-4e1d-a6d3-0ed64a449490
# ╠═67c63d90-496a-4fca-9146-b2ba24ebb4be
