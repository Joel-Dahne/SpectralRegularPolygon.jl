### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ cc7c7dea-83ab-11f0-04b8-a7eb5cc4403a
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

    # We use the SRP.S function a lot, so import this separately
    using SpectralRegularPolygons: S

    setprecision(Arb, 128)
end

# ╔═╡ c45a43db-337d-4ce7-afeb-61a9f50ba398
md"""
# Proof of Lemma 2.6
"""

# ╔═╡ 4de76364-83fd-4889-91ff-cf0ae92ffdea
md"""
## Goal
We want to prove that for $|z| = 1$ we have the following bounds:

$$|b_2(z)| \leq C_{b,2},\ |b_3(z)| \leq C_{b,3},\ |b_4(z)| \leq C_{b,4},\ |b_5(z)| \leq C_{b,5},$$

and for $N \geq N_0$ have

$$|T_2(z)| \leq C_{T,2},\ |T_4(z)| \leq C_{T,4},\ |T_6(z)| \leq C_{T,6}.$$

Here $N_0$, $C_{b,k}$ and $C_{T,l}$ are given by:
"""

# ╔═╡ b20ebf73-ffca-41f6-a2be-13d56881ad25
N₀ = SRP.N₀

# ╔═╡ bff31288-f074-4e51-a4b5-b97de10af009
C_b_2 = SRP.C_b_2

# ╔═╡ 53e1c0b9-ae24-49a9-90a9-a7d2c9aa0a34
C_b_3 = SRP.C_b_3

# ╔═╡ cb12a4fa-c971-4551-8481-dfd35c1f824d
C_b_4 = SRP.C_b_4

# ╔═╡ 0840b31c-ae3c-408f-a991-16828560178c
C_b_5 = SRP.C_b_5

# ╔═╡ ea0e8a40-c152-4be3-b442-cfc85fb8f3e9
C_T_2 = SRP.C_T_2

# ╔═╡ cba219da-6ae9-465b-8f47-7681945c8421
C_T_4 = SRP.C_T_4

# ╔═╡ ae26b910-dc30-4904-aa25-c9f0a860d6ca
C_T_6 = SRP.C_T_6

# ╔═╡ 9a1e011a-6875-4a11-b321-b703a2e15e62
md"""
## Implementations of functions
The functions $b_2$, $b_3$, $b_4$ and $b_5$ are straight forward to implement:
"""

# ╔═╡ f730437f-c242-4036-9688-85efe213a398
b_2(z) = real(S(2, z))

# ╔═╡ b4167ebb-b588-45fe-acc4-0653ed6d2a9a
b_3(z) = real(S(3, z))

# ╔═╡ 88e1d193-ef1d-44de-a9ca-a8bc66a5dd0f
b_4(z) = imag(S(2, z))^2 / 2 + real(S(4, z))

# ╔═╡ f9ea139f-bf69-4919-8192-460c668b1bf3
b_5(z) = imag(S(2, z)) * imag(S(3, z)) + real(S(5, z)) - SRP.λ_disc() * zeta(Arb(5))

# ╔═╡ 0097219d-cfb9-4aa2-a6d1-6f1d232dc7eb
md"""
To compute $T_6$ we first compute a Taylor model for

$$\frac{\sqrt{\lambda_{\text{app}}}}{\sqrt{\lambda}} c_N |F_N(z)|.$$

in terms of $N^{-1}$. If we compute the Taylor model with a remainder term of degree 6, then $T_6$ is bounded by this remainder term.

Since the factor

$$\frac{\sqrt{\lambda_{\text{app}}}}{\sqrt{\lambda}} c_N$$

doesn't depend on $z$ we can precompute it:
"""

# ╔═╡ fa92b112-5b49-456c-b207-8803a43a10fd
# It is really important to get a good enclosure of this factor,
# so we compute a higher order Taylor model and then truncate it.
λ_app_div_λ_mul_c_N = SRP.truncate(
    SRP.ArbTaylorModel(
        inv_N -> sqrt(SRP.λ_app_div_λ(inv_N)) * SRP.c_N(inv_N),
        Arb((0, 1 // N₀)),
        Arb(0),
        degree = 10,
    ),
    degree = 5,
)

# ╔═╡ b4186719-9e4b-4d24-9c26-adcb33acc3dc
md"""
TODO: Explain this a bit more.
"""

# ╔═╡ 26b77fc2-2e84-4d8f-8cbc-6a15a948aab7
T_6(z) = (λ_app_div_λ_mul_c_N * abs(SRP.F_N_model(N₀, z))).p[end]

# ╔═╡ e707c5eb-acd8-43fc-8fd1-1e0f9683d96b
md"""
## Plots
"""

# ╔═╡ 0aa4c803-768e-4b65-ba2c-c17f400d17d9
θs_div_πs = range(Arb(0), 2, 1000)

# ╔═╡ 28ef9722-2d9c-4ad3-ac78-baec417e9a45
θs = π * θs_div_πs

# ╔═╡ 5d2dd24d-e9bd-434a-9a34-c4b9c130f1e6
b2s = tmap(θs_div_πs) do θ_div_π
    b_2(SRP.exppii(θ_div_π))
end

# ╔═╡ e0ccdc24-1cd5-4967-86a4-226423d33fb5
b3s = tmap(θs_div_πs) do θ_div_π
    b_3(SRP.exppii(θ_div_π))
end

# ╔═╡ f8466db6-f0d4-44cc-b369-a587f8565bc7
b4s = tmap(θs_div_πs) do θ_div_π
    b_4(SRP.exppii(θ_div_π))
end

# ╔═╡ 95b752e2-5927-4235-babb-5a2198f172aa
b5s = tmap(θs_div_πs) do θ_div_π
    b_5(SRP.exppii(θ_div_π))
end

# ╔═╡ c25aea0b-f364-4c25-9fa8-8315fcc319cb
T6s = tmap(θs_div_πs) do θ_div_π
    T_6(SRP.exppii(θ_div_π))
end

# ╔═╡ 0209a6b1-f5ba-4ba6-9b5f-2f2cbdd2c643
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"b_2(e^{i\theta})")
    band!(ax, θs, lbound.(b2s), ubound.(b2s))
    scatterlines!(ax, θs, b2s)
    hlines!(ax, [-Arb(C_b_2), Arb(C_b_2)])
    fig
end

# ╔═╡ 1441945f-d537-4ac6-b21d-cfcc62c28cef
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"b_3(e^{i\theta})")
    band!(ax, θs, lbound.(b3s), ubound.(b3s))
    scatterlines!(ax, θs, b3s)
    hlines!(ax, [-Arb(C_b_3), Arb(C_b_3)])
    fig
end

# ╔═╡ 0c002f46-d981-4461-a2fe-0ceaf79ece35
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"b_4(e^{i\theta})")
    band!(ax, θs, lbound.(b4s), ubound.(b4s))
    scatterlines!(ax, θs, b4s)
    hlines!(ax, [-Arb(C_b_4), Arb(C_b_4)])
    fig
end

# ╔═╡ 0d186cc4-b8cd-4336-aa8d-e2ffdc37c80a
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"b_5(e^{i\theta})")
    band!(ax, θs, lbound.(b5s), ubound.(b5s))
    scatterlines!(ax, θs, b5s)
    hlines!(ax, [-Arb(C_b_5), Arb(C_b_5)])
    fig
end

# ╔═╡ 1f4ebea3-1b83-49d4-aabc-cc57b0bba1af
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\theta", ylabel = L"T_6(e^{i\theta})")
    band!(ax, θs, lbound.(T6s), ubound.(T6s))
    scatterlines!(ax, θs, T6s)
    hlines!(ax, [-Arb(C_T_6), Arb(C_T_6)])
    fig
end

# ╔═╡ e6f37ebc-4884-4028-a937-93d71c27e3d0
md"""
## Proofs
"""

# ╔═╡ 0a8e3984-1047-4f89-8f16-d47b89fcd6cf
@time b_2_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    b_2(SRP.exppii(θ_div_π))
end

# ╔═╡ fd627d82-96e2-4f57-923b-dd879acf5bc6
@time b_3_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    b_3(SRP.exppii(θ_div_π))
end

# ╔═╡ b3f8fb15-2864-4479-8b2a-74586b2c8e8b
@time b_4_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    b_4(SRP.exppii(θ_div_π))
end

# ╔═╡ 00ac7282-9dac-4289-8d5a-5bf5dd539878
@time b_5_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    b_5(SRP.exppii(θ_div_π))
end

# ╔═╡ 763878e0-5a56-4c0e-985f-f0521cab30eb
@time T_6_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    ubound_tol = Arb(C_T_6),
    depth = 30,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    T_6(SRP.exppii(θ_div_π))
end

# ╔═╡ 92e50091-4405-4691-aec5-8704fac4215a
getinterval(Arb, T_6_bound) # TODO: Format this in a nice way

# ╔═╡ dad6497c-bfee-4235-82e6-2099220844b3
T_2_bound =
    b_2_bound + b_3_bound / N₀ + b_4_bound / N₀^2 + b_5_bound / N₀^3 + T_6_bound / N₀^4

# ╔═╡ 769879b3-3c7b-4ae5-ad75-5f5fee225f90
T_4_bound = b_4_bound + b_5_bound / N₀ + T_6_bound / N₀^2

# ╔═╡ 48596d0e-6804-426a-b3a4-7d64f3d7df93
b_2_bound <= Arb(C_b_2)

# ╔═╡ b9322190-84d9-4a11-b708-8ff8c11fc367
b_3_bound <= Arb(C_b_3)

# ╔═╡ 7ad99c3c-5dfe-4518-bf91-113388af8a2c
b_4_bound <= Arb(C_b_4)

# ╔═╡ 8f495d5a-03cd-4d31-ad19-2fb4fde08d70
b_5_bound <= Arb(C_b_5)

# ╔═╡ a9cb0e9a-f69e-45d0-995a-cea5db9c22ec
T_2_bound <= Arb(C_T_2)

# ╔═╡ d9af7e04-5297-4a4f-961c-171b2c939d8f
T_4_bound <= Arb(C_T_4)

# ╔═╡ 26926dfc-0236-4f86-a7be-a8881e6f5f30
T_6_bound <= Arb(C_T_6)

# ╔═╡ Cell order:
# ╟─c45a43db-337d-4ce7-afeb-61a9f50ba398
# ╠═cc7c7dea-83ab-11f0-04b8-a7eb5cc4403a
# ╟─4de76364-83fd-4889-91ff-cf0ae92ffdea
# ╠═b20ebf73-ffca-41f6-a2be-13d56881ad25
# ╠═bff31288-f074-4e51-a4b5-b97de10af009
# ╠═53e1c0b9-ae24-49a9-90a9-a7d2c9aa0a34
# ╠═cb12a4fa-c971-4551-8481-dfd35c1f824d
# ╠═0840b31c-ae3c-408f-a991-16828560178c
# ╠═ea0e8a40-c152-4be3-b442-cfc85fb8f3e9
# ╠═cba219da-6ae9-465b-8f47-7681945c8421
# ╠═ae26b910-dc30-4904-aa25-c9f0a860d6ca
# ╟─9a1e011a-6875-4a11-b321-b703a2e15e62
# ╠═f730437f-c242-4036-9688-85efe213a398
# ╠═b4167ebb-b588-45fe-acc4-0653ed6d2a9a
# ╠═88e1d193-ef1d-44de-a9ca-a8bc66a5dd0f
# ╠═f9ea139f-bf69-4919-8192-460c668b1bf3
# ╟─0097219d-cfb9-4aa2-a6d1-6f1d232dc7eb
# ╠═fa92b112-5b49-456c-b207-8803a43a10fd
# ╟─b4186719-9e4b-4d24-9c26-adcb33acc3dc
# ╠═26b77fc2-2e84-4d8f-8cbc-6a15a948aab7
# ╟─e707c5eb-acd8-43fc-8fd1-1e0f9683d96b
# ╠═0aa4c803-768e-4b65-ba2c-c17f400d17d9
# ╠═28ef9722-2d9c-4ad3-ac78-baec417e9a45
# ╠═5d2dd24d-e9bd-434a-9a34-c4b9c130f1e6
# ╠═e0ccdc24-1cd5-4967-86a4-226423d33fb5
# ╠═f8466db6-f0d4-44cc-b369-a587f8565bc7
# ╠═95b752e2-5927-4235-babb-5a2198f172aa
# ╠═c25aea0b-f364-4c25-9fa8-8315fcc319cb
# ╟─0209a6b1-f5ba-4ba6-9b5f-2f2cbdd2c643
# ╟─1441945f-d537-4ac6-b21d-cfcc62c28cef
# ╟─0c002f46-d981-4461-a2fe-0ceaf79ece35
# ╟─0d186cc4-b8cd-4336-aa8d-e2ffdc37c80a
# ╠═1f4ebea3-1b83-49d4-aabc-cc57b0bba1af
# ╟─e6f37ebc-4884-4028-a937-93d71c27e3d0
# ╠═0a8e3984-1047-4f89-8f16-d47b89fcd6cf
# ╠═fd627d82-96e2-4f57-923b-dd879acf5bc6
# ╠═b3f8fb15-2864-4479-8b2a-74586b2c8e8b
# ╠═00ac7282-9dac-4289-8d5a-5bf5dd539878
# ╠═763878e0-5a56-4c0e-985f-f0521cab30eb
# ╠═92e50091-4405-4691-aec5-8704fac4215a
# ╠═dad6497c-bfee-4235-82e6-2099220844b3
# ╠═769879b3-3c7b-4ae5-ad75-5f5fee225f90
# ╠═48596d0e-6804-426a-b3a4-7d64f3d7df93
# ╠═b9322190-84d9-4a11-b708-8ff8c11fc367
# ╠═7ad99c3c-5dfe-4518-bf91-113388af8a2c
# ╠═8f495d5a-03cd-4d31-ad19-2fb4fde08d70
# ╠═a9cb0e9a-f69e-45d0-995a-cea5db9c22ec
# ╠═d9af7e04-5297-4a4f-961c-171b2c939d8f
# ╠═26926dfc-0236-4f86-a7be-a8881e6f5f30
