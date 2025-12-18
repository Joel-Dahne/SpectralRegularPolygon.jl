### A Pluto.jl notebook ###
# v0.20.21

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

    setprecision(Arb, 128)
end

# ╔═╡ 18d792d6-8c0c-4214-bd98-670611da3e1e
md"""
# Proof of Lemma 2.10
"""

# ╔═╡ bdc6dcff-50a0-49d7-98a2-00b24f5ecfe6
md"""
We want to prove that for $|z| = 1$ we have the following bounds:

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} d_k(z, t)V_l(t)\ dt\right| \leq C_{I,k,l}$$

for

$$(k, l) = (1, 4), (2, 3), (2, 4), (3, 2), (3, 3) \text{ and } (3, 4).$$

For $N \geq N_0$ we want to prove that we have

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} K_4(z, t)V_l(t)\ dt\right| \leq C_{I,K,l}$$

for $l = 1, 2, 3, 4$.

For the cases when $k + l > 5$ in the first case and $l > 1$ in the second case we we only need very rough bounds. These constants are eventually divided by some power of $N_0$ and therefore do not to be as precise. 

Here $N_0$$, $C_{I,k,l}$ and $C_{K,l}$ are given by
"""

# ╔═╡ d22dc5aa-bac2-4f8a-968b-ece39e295129
N₀ = SRP.N₀

# ╔═╡ 5c5f73da-4135-412f-a31f-da02450f235a
C_I_1_4 = SRP.C_I_1_4

# ╔═╡ 4b108a4d-cd35-4f5c-ad2f-523d6c3c73d9
C_I_2_3 = SRP.C_I_2_3

# ╔═╡ 1d9f8eb1-e9b1-4ffa-9c7d-19567d582e85
C_I_2_4 = SRP.C_I_2_4

# ╔═╡ 11392f25-b796-4f5d-a970-c0a2f82bdcd6
C_I_3_2 = SRP.C_I_3_2

# ╔═╡ 978b5167-ca3f-4717-965f-e71e55280997
C_I_3_3 = SRP.C_I_3_3

# ╔═╡ 1c8ee061-1d61-4e58-a206-20466ff20b63
C_I_3_4 = SRP.C_I_3_4

# ╔═╡ 278f9795-6b9b-4780-a0c7-87a12efa8e09
C_I_K_1 = SRP.C_I_K_1

# ╔═╡ 37b7d3ac-e62e-4dbb-8790-bbc285dc710b
C_I_K_2 = SRP.C_I_K_2

# ╔═╡ c089ba08-10c7-4df2-9b0d-cf1be9b9d2ea
C_I_K_3 = SRP.C_I_K_3

# ╔═╡ 6c4607e4-b04e-4698-adee-446d331c4594
C_I_K_4 = SRP.C_I_K_4

# ╔═╡ b6096ef5-2ce3-4519-97db-6ccb566037f0
md"""
## Plots
To get a better understanding for how the functions behave we plot them on the interval $[0, 2\pi]$. Note that these plots are not part of the proof, they are only mean to give an idea for what the functions behave like.

For the case when $k + l > 5$ for the integrals with $d$ or $l > 1$ for the integrals with $K_4$ we compute very rough enclosures and the plots below do not give much information about the actual behavior of the function. They do however give an indication of the type of bounds we could prove.
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

# ╔═╡ 5844760d-3dbe-42bd-b5bd-e49fd94ea01e
I_d_2_V_4 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_d_k_V_l(2, 4, SRP.exppii(θ_div_π))
end

# ╔═╡ a820d24b-9387-4ef5-a959-a4baaf93fdb2
I_d_3_V_2 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_d_k_V_l(3, 2, SRP.exppii(θ_div_π))
end

# ╔═╡ 952ad264-430b-4ab4-9c47-c6af0bbb1de6
I_d_3_V_3 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_d_k_V_l(3, 3, SRP.exppii(θ_div_π))
end

# ╔═╡ 5f8946a7-3da8-4786-9248-28c90f482db3
I_d_3_V_4 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_d_k_V_l(3, 4, SRP.exppii(θ_div_π))
end

# ╔═╡ a0bcac27-bbf7-481e-96bd-f1b619e7a54a
I_K_4_V_1 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_K_4_V_l(N₀, 1, SRP.exppii(θ_div_π))
end

# ╔═╡ 92c7561b-b17f-4138-826b-eef30c8226f8
I_K_4_V_2 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_K_4_V_l(N₀, 2, SRP.exppii(θ_div_π))
end

# ╔═╡ 4b0a4d2c-c23c-416d-9356-108f1ab1d17b
I_K_4_V_3 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_K_4_V_l(N₀, 3, SRP.exppii(θ_div_π))
end

# ╔═╡ 12c3033e-5b6c-4967-b339-dad84f4806e9
I_K_4_V_4 = tmap(θs_div_πs) do θ_div_π
    SRP.integral_K_4_V_l(N₀, 4, SRP.exppii(θ_div_π))
end

# ╔═╡ d102ccd5-5db0-4d3c-93b3-7e3b6da4c958
md"""
In the below plots, the horizontal lines indicate the bounds that are proved to hold.
"""

# ╔═╡ 7bc1266a-4f01-439b-81fb-5b1e86612ef2
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{1,4}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_d_1_V_4)), ubound.(real(I_d_1_V_4)))
    hlines!(ax, [-Arb(C_I_1_4), Arb(C_I_1_4)])
    fig
end

# ╔═╡ 6187163b-832d-4c9f-862f-a3483505fda6
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{2,3}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_d_2_V_3)), ubound.(real(I_d_2_V_3)))
    hlines!(ax, [-Arb(C_I_2_3), Arb(C_I_2_3)])
    fig
end

# ╔═╡ 53298c0c-b59d-4e38-8c27-ef1f755f6e51
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{2,4}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_d_2_V_4)), ubound.(real(I_d_2_V_4)))
    hlines!(ax, [-Arb(C_I_2_4), Arb(C_I_2_4)])
    fig
end

# ╔═╡ a7509228-24ea-4395-8bdc-57457c63e28b
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{3,2}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_d_3_V_2)), ubound.(real(I_d_3_V_2)))
    hlines!(ax, [-Arb(C_I_3_2), Arb(C_I_3_2)])
    fig
end

# ╔═╡ b18bac31-c245-4d8e-9bbe-11f85678714c
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{3,3}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_d_3_V_3)), ubound.(real(I_d_3_V_3)))
    hlines!(ax, [-Arb(C_I_3_3), Arb(C_I_3_3)])
    fig
end

# ╔═╡ d68d454f-baf3-4b36-914f-f3a81226da69
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{3,4}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_d_3_V_4)), ubound.(real(I_d_3_V_4)))
    hlines!(ax, [-Arb(C_I_3_4), Arb(C_I_3_4)])
    fig
end

# ╔═╡ bab7bd70-2f53-4421-be57-c79c469c6a5a
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{K,1}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_K_4_V_1)), ubound.(real(I_K_4_V_1)))
    hlines!(ax, [-Arb(C_I_K_1), Arb(C_I_K_1)])
    fig
end

# ╔═╡ a5a011f7-f447-4ad6-a0f8-ee6a18490cf2
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{K,2}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_K_4_V_2)), ubound.(real(I_K_4_V_2)))
    hlines!(ax, [-Arb(C_I_K_2), Arb(C_I_K_2)])
    fig
end

# ╔═╡ 93398ee3-2183-4f26-937b-ee40cac51db7
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{K,3}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_K_4_V_3)), ubound.(real(I_K_4_V_3)))
    hlines!(ax, [-Arb(C_I_K_3), Arb(C_I_K_3)])
    fig
end

# ╔═╡ cd7af643-b584-4f24-b058-05738aea9284
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"I_{K,4}(e^{i\theta})")
    band!(ax, θs, lbound.(real(I_K_4_V_4)), ubound.(real(I_K_4_V_4)))
    hlines!(ax, [-Arb(C_I_K_4), Arb(C_I_K_4)])
    fig
end

# ╔═╡ 93cf97f7-54c1-4101-93b8-09890db8de86
md"""
## Proof
"""

# ╔═╡ f293357f-b59f-46e8-9b68-a83f1a7cbe42
#=╠═╡
@time I_d_1_V_4_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    ubound_tol = Arb(C_I_1_4),
    depth_start = 8,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_d_k_V_l(1, 4, SRP.exppii(θ_div_π)))
end
  ╠═╡ =#

# ╔═╡ 0cfdcd21-461a-4626-b7c8-f6257b9ff3e6
@time I_d_2_V_3_bound = ArbExtras.maximum_enclosure(
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

# ╔═╡ ac42d8a6-2373-40bb-acd9-8010caabbe14
@time I_d_2_V_4_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    ubound_tol = Arb(C_I_2_4),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_d_k_V_l(2, 4, SRP.exppii(θ_div_π)))
end

# ╔═╡ 3d283f17-149d-4e1d-a6d3-0ed64a449490
@time I_d_3_V_2_bound = ArbExtras.maximum_enclosure(
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

# ╔═╡ c018a526-d2ac-405b-a3ac-5b267095cb52
@time I_d_3_V_3_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    ubound_tol = Arb(C_I_3_3),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_d_k_V_l(3, 3, SRP.exppii(θ_div_π)))
end

# ╔═╡ ea0b645a-c745-459a-b29d-8f76b9edf3ac
@time I_d_3_V_4_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    ubound_tol = Arb(C_I_3_4),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_d_k_V_l(3, 4, SRP.exppii(θ_div_π)))
end

# ╔═╡ 67c63d90-496a-4fca-9146-b2ba24ebb4be
@time I_K_4_V_1_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    ubound_tol = Arb(C_I_K_1),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_K_4_V_l(N₀, 1, SRP.exppii(θ_div_π)))
end

# ╔═╡ 583a4acf-52e9-4aba-9cca-2e12865ca30e
@time I_K_4_V_2_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    ubound_tol = Arb(C_I_K_2),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_K_4_V_l(N₀, 2, SRP.exppii(θ_div_π)))
end

# ╔═╡ 849fe7d7-bc52-48f0-96e0-590332f76694
@time I_K_4_V_3_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    ubound_tol = Arb(C_I_K_3),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_K_4_V_l(N₀, 3, SRP.exppii(θ_div_π)))
end

# ╔═╡ 61037aaa-08d2-46eb-b206-b89f94ef48fd
@time I_K_4_V_4_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    ubound_tol = Arb(C_I_K_4),
    depth_start = 4,
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    real(SRP.integral_K_4_V_l(N₀, 4, SRP.exppii(θ_div_π)))
end

# ╔═╡ 2ff324b1-da1d-446f-89f5-580aad93681d
md"""
For $I_{2,3}$ and $I_K$ the bounds are faily wide and printing them as balls is not optimal. We therefore print them as intervals with lower and upper bounds.
"""

# ╔═╡ 4b598e95-0b9b-4d45-863f-b4841959664c
string.(getinterval(Arb, I_d_2_V_3_bound), digits = 10)

# ╔═╡ f632c4a3-6216-49e5-a49d-0c92e9048309
string.(getinterval(Arb, I_K_4_V_1_bound), digits = 10)

# ╔═╡ c73dbe22-cacf-42c7-96b0-017682087acc
md"""
Finally we verify that the required bounds are satisfied.
"""

# ╔═╡ 2588eeda-711f-4acd-91bf-24349dcd90ba
@assert_proof I_d_1_V_4_bound <= Arb(C_I_1_4)

# ╔═╡ 9d1c43d5-3c18-4c78-94c9-df9f6e0912c5
@assert_proof I_d_2_V_3_bound <= Arb(C_I_2_3)

# ╔═╡ 41e47cbf-1c0f-4e78-a375-7f53aaa9af03
@assert_proof I_d_2_V_4_bound <= Arb(C_I_2_4)

# ╔═╡ 73f4c9c7-0975-4a0a-a0cc-673f449ad8c7
@assert_proof I_d_3_V_2_bound <= Arb(C_I_3_2)

# ╔═╡ 3a232ec8-68a6-4aab-b209-e3d1aa58061a
@assert_proof I_d_3_V_3_bound <= Arb(C_I_3_3)

# ╔═╡ 73bd865a-b043-46f5-bc21-812bd16cb2b6
@assert_proof I_d_3_V_4_bound <= Arb(C_I_3_4)

# ╔═╡ b11334b3-19f8-4dbc-bb70-dd137d342104
@assert_proof I_K_4_V_1_bound <= Arb(C_I_K_1)

# ╔═╡ 4a30bf9e-1e4e-4f69-b8ec-5ab3844b1023
@assert_proof I_K_4_V_2_bound <= Arb(C_I_K_2)

# ╔═╡ 7ceec59e-a3ff-4ef3-9593-0555ff6dc6b3
@assert_proof I_K_4_V_3_bound <= Arb(C_I_K_3)

# ╔═╡ 7d1720c5-016c-44bc-b7d4-ea2f3a3eed48
@assert_proof I_K_4_V_4_bound <= Arb(C_I_K_4)

# ╔═╡ Cell order:
# ╟─18d792d6-8c0c-4214-bd98-670611da3e1e
# ╠═fca89648-853b-11f0-3116-5b610b900eb4
# ╟─bdc6dcff-50a0-49d7-98a2-00b24f5ecfe6
# ╠═d22dc5aa-bac2-4f8a-968b-ece39e295129
# ╠═5c5f73da-4135-412f-a31f-da02450f235a
# ╠═4b108a4d-cd35-4f5c-ad2f-523d6c3c73d9
# ╠═1d9f8eb1-e9b1-4ffa-9c7d-19567d582e85
# ╠═11392f25-b796-4f5d-a970-c0a2f82bdcd6
# ╠═978b5167-ca3f-4717-965f-e71e55280997
# ╠═1c8ee061-1d61-4e58-a206-20466ff20b63
# ╠═278f9795-6b9b-4780-a0c7-87a12efa8e09
# ╠═37b7d3ac-e62e-4dbb-8790-bbc285dc710b
# ╠═c089ba08-10c7-4df2-9b0d-cf1be9b9d2ea
# ╠═6c4607e4-b04e-4698-adee-446d331c4594
# ╟─b6096ef5-2ce3-4519-97db-6ccb566037f0
# ╠═ef2ecb67-904b-4488-8e47-b5e2fc46eaab
# ╠═d24e7344-a2cd-4ba3-a5cf-1290fe285cfc
# ╠═cf4e04c7-26f2-400b-a885-c3e495e42f78
# ╠═67e1ad91-7da5-4985-b8c3-ca1cc84a47fb
# ╠═5844760d-3dbe-42bd-b5bd-e49fd94ea01e
# ╠═a820d24b-9387-4ef5-a959-a4baaf93fdb2
# ╠═952ad264-430b-4ab4-9c47-c6af0bbb1de6
# ╠═5f8946a7-3da8-4786-9248-28c90f482db3
# ╠═a0bcac27-bbf7-481e-96bd-f1b619e7a54a
# ╠═92c7561b-b17f-4138-826b-eef30c8226f8
# ╠═4b0a4d2c-c23c-416d-9356-108f1ab1d17b
# ╠═12c3033e-5b6c-4967-b339-dad84f4806e9
# ╟─d102ccd5-5db0-4d3c-93b3-7e3b6da4c958
# ╟─7bc1266a-4f01-439b-81fb-5b1e86612ef2
# ╟─6187163b-832d-4c9f-862f-a3483505fda6
# ╟─53298c0c-b59d-4e38-8c27-ef1f755f6e51
# ╟─a7509228-24ea-4395-8bdc-57457c63e28b
# ╟─b18bac31-c245-4d8e-9bbe-11f85678714c
# ╟─d68d454f-baf3-4b36-914f-f3a81226da69
# ╟─bab7bd70-2f53-4421-be57-c79c469c6a5a
# ╟─a5a011f7-f447-4ad6-a0f8-ee6a18490cf2
# ╟─93398ee3-2183-4f26-937b-ee40cac51db7
# ╟─cd7af643-b584-4f24-b058-05738aea9284
# ╟─93cf97f7-54c1-4101-93b8-09890db8de86
# ╠═f293357f-b59f-46e8-9b68-a83f1a7cbe42
# ╠═0cfdcd21-461a-4626-b7c8-f6257b9ff3e6
# ╠═ac42d8a6-2373-40bb-acd9-8010caabbe14
# ╠═3d283f17-149d-4e1d-a6d3-0ed64a449490
# ╠═c018a526-d2ac-405b-a3ac-5b267095cb52
# ╠═ea0b645a-c745-459a-b29d-8f76b9edf3ac
# ╠═67c63d90-496a-4fca-9146-b2ba24ebb4be
# ╠═583a4acf-52e9-4aba-9cca-2e12865ca30e
# ╠═849fe7d7-bc52-48f0-96e0-590332f76694
# ╠═61037aaa-08d2-46eb-b206-b89f94ef48fd
# ╟─2ff324b1-da1d-446f-89f5-580aad93681d
# ╠═4b598e95-0b9b-4d45-863f-b4841959664c
# ╠═f632c4a3-6216-49e5-a49d-0c92e9048309
# ╟─c73dbe22-cacf-42c7-96b0-017682087acc
# ╠═2588eeda-711f-4acd-91bf-24349dcd90ba
# ╠═9d1c43d5-3c18-4c78-94c9-df9f6e0912c5
# ╠═41e47cbf-1c0f-4e78-a375-7f53aaa9af03
# ╠═73f4c9c7-0975-4a0a-a0cc-673f449ad8c7
# ╠═3a232ec8-68a6-4aab-b209-e3d1aa58061a
# ╠═73bd865a-b043-46f5-bc21-812bd16cb2b6
# ╠═b11334b3-19f8-4dbc-bb70-dd137d342104
# ╠═4a30bf9e-1e4e-4f69-b8ec-5ab3844b1023
# ╠═7ceec59e-a3ff-4ef3-9593-0555ff6dc6b3
# ╠═7d1720c5-016c-44bc-b7d4-ea2f3a3eed48
