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
    using PlutoUI
    using SpecialFunctions

    import SpectralRegularPolygons as SRP

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

$$\frac{\lambda_{\text{app}}(N)}{1 + \hat{\varepsilon}(N)} > \frac{\lambda_{\text{app}}(N+1)}{1 - \hat{\varepsilon}(N+1)}$$

and

$$\frac{\lambda_{\text{app}}(N)(1-\hat{\varepsilon}(N+1))}{(1+\hat{\varepsilon}(N))\lambda_{\text{app}}(N+1)} > \frac{\lambda_{\text{app}}(N+1)(1+\hat{\varepsilon}(N+2))}{(1-\hat{\varepsilon}(N+1))\lambda_{\text{app}}(N+2)},$$

where $N_0$ is given by:
"""

# ╔═╡ 1e814a20-d22d-4c0d-9003-6712a0a5d205
N₀ = SRP.N₀

# ╔═╡ ecdfb51c-bd88-4d2b-86fd-d6c46340bf28
md"""
## Proof
### Part 1
Let us start by proving that

$$\frac{\lambda_{\text{app}}(N)}{1 + \hat{\varepsilon}(N)} > \frac{\lambda_{\text{app}}(N+1)}{1 - \hat{\varepsilon}(N+1)}.$$

It is convenient to see the inequality as a function of $N^{-1}$, rather than $N$. For that reason we introduce the function

$$G(N^{-1}) = \frac{\lambda_{\text{app}}(N)}{1 + \hat{\varepsilon}(N)} - \frac{\lambda_{\text{app}}(N+1)}{1 - \hat{\varepsilon}(N+1)}.$$
"""

# ╔═╡ 03f4726b-2d20-4a77-83ba-41824798533f
function G(inv_N)
    inv_Np1 = inv_N / (1 + inv_N) # Enclosure of inv(N + 1)

    return SRP.λ_disc() * (
        SRP.λ_app_div_λ(inv_N) / (1 + SRP.epsilon_hat(inv_N)) -
        SRP.λ_app_div_λ(inv_Np1) / (1 - SRP.epsilon_hat(inv_Np1))
    )
end

# ╔═╡ bb882554-bd78-4ab2-acfc-e7fb702e0024
md"""
The goal is then to prove that $G(N^{-1}) > 0$ for $0 < N^{-1} \leq N_0^{-1}$.

Note that $G(0) = 0$, so we cannot directly enclose $G$ and prove that it is zero. Furthermore, $G'(0) = G''(0) = G'''(0) = 0$ but $G''''(0) > 0$. We can see this by computing the Taylor expansion of $G$ at zero:
"""

# ╔═╡ 76af3933-de21-45f0-9310-028eaf31465f
G_expansion = G(ArbSeries((0, 1), degree = 4))

# ╔═╡ 4eb02435-b201-4e16-9ac5-fbc68c131dad
md"""
We immediately get that $G(0) = G'(0) = G''(0) = 0$:
"""

# ╔═╡ b20f3e75-b475-4cac-b69d-b6611ac82363
iszero(G_expansion[0]) && iszero(G_expansion[1]) && iszero(G_expansion[2])

# ╔═╡ bfcf8109-a214-4460-b41a-35a714096f05
md"""
For $G'''(0)$ we get an enclosure containing zero:
"""

# ╔═╡ 0d216fba-d53a-47ba-886f-dc8902b0b0c2
Arblib.contains_zero(G_expansion[3])

# ╔═╡ e03b99aa-2802-470e-adb3-dbf6a74da074
md"""
But this doesn't prove that it is zero, this is proved in the paper.

Since $G'(0) = G''(0) = G'''(0) = 0$, it suffices to prove that $G^{(4)}(N^{-1}) > 0$ for $0 \leq N^{-1} \leq N_0^{-1}$. We therefore enclose the minimum value of $G^{(4)}(N^{-1})$ on this interval:
"""

# ╔═╡ f7c11696-e359-48c6-bf29-b799b3a26957
G_d4_minimum = ArbExtras.minimum_enclosure(
    ArbExtras.derivative_function(G, 4),
    Arf(0),
    ubound(Arb(1 // N₀)),
    verbose = true,
)

# ╔═╡ b2b04733-b22f-4283-b657-013bd02864ac
md"""
Finally we verify that the minimum is positive:
"""

# ╔═╡ 1b5615d5-ff94-41c6-8123-d11bdae8af08
Arblib.ispositive(G_d4_minimum)

# ╔═╡ 8184b4d2-7235-4cfe-8a34-a267ae170be7
md"""
### Part 2
"""

# ╔═╡ 4371e67c-6c2e-461c-bbb5-afcb8ad7f142
md"""
Next we prove that

$$\frac{\lambda_{\text{app}}(N)(1-\hat{\varepsilon}(N+1))}{(1+\hat{\varepsilon}(N))\lambda_{\text{app}}(N+1)} > \frac{\lambda_{\text{app}}(N+1)(1+\hat{\varepsilon}(N+2))}{(1-\hat{\varepsilon}(N+1))\lambda_{\text{app}}(N+2)}.$$

The approach is similar to the one above, we let

$$G_q(N^{-1}) = \frac{\lambda_{\text{app}}(N)(1-\hat{\varepsilon}(N+1))}{(1+\hat{\varepsilon}(N))\lambda_{\text{app}}(N+1)} - \frac{\lambda_{\text{app}}(N+1)(1+\hat{\varepsilon}(N+2))}{(1-\hat{\varepsilon}(N+1))\lambda_{\text{app}}(N+2)}.$$
"""

# ╔═╡ 89947765-e753-4750-b6cc-084c0c86280f
function G_q(inv_N)
    # Enclosures of inv(N + 1) and inv(N + 2)
    inv_Np1 = inv_N / (1 + inv_N)
    inv_Np2 = inv_Np1 / (1 + inv_Np1)

    # Note that we use λ_app_div_λ since the λs cancel
    return SRP.λ_app_div_λ(inv_N) * (1 - SRP.epsilon_hat(inv_Np1)) /
           ((1 + SRP.epsilon_hat(inv_N)) * SRP.λ_app_div_λ(inv_Np1)) -
           SRP.λ_app_div_λ(inv_Np1) * (1 + SRP.epsilon_hat(inv_Np2)) /
           ((1 - SRP.epsilon_hat(inv_Np1)) * SRP.λ_app_div_λ(inv_Np2))
end

# ╔═╡ 80973d1e-71af-454e-b0ee-a949e445ed44
md"""
The goal is then to prove that $G_q(N^{-1}) > 0$ for $0 < N^{-1} \leq N_0^{-1}$.

Similar to for $G$ we have $G_q(0) = G_q'(0) = G_q''(0) = G_q'''(0) = 0$ and in this case we also have $G_q^{(4)}(0) = 0$, but $G_q^{(5)}(0) > 0$. We can see this by computing the Taylor expansion of $G_q$ at zero:
"""

# ╔═╡ c427a061-0bf9-4e40-947e-2bd093ad686a
G_q_expansion = G_q(ArbSeries((0, 1), degree = 5))

# ╔═╡ 40f8bc3e-9d69-4e95-800c-7f33e5187b8c
md"""
We immediately get that $G_q(0) = G_q'(0) = G_q''(0) = 0$:
"""

# ╔═╡ 8374aeda-b25a-4a9a-b392-47b423892b86
iszero(G_q_expansion[0]) && iszero(G_q_expansion[1]) && iszero(G_q_expansion[2])

# ╔═╡ 844d9803-36ca-40a2-8dbf-911a52ede8dc
md"""
For $G_q'''(0)$ and $G_q^{(4)}(0)$ we get enclosures containing zero:
"""

# ╔═╡ 218d06e7-c755-41b3-bfaa-1dfd6a3ccded
Arblib.contains_zero(G_q_expansion[3]) && Arblib.contains_zero(G_q_expansion[4])

# ╔═╡ bcbcb115-b950-4e71-ab50-665e67b32c68
md"""
But this doesn't prove they are zero, this is proved in the paper.

Contrary to for $G$ we do not have $G_q^{(5)}(N^{-1}) > 0$ for the entire interval $0 \leq N^{-1} \leq N_0^{-1}$, but only part of it. We therefore split the interval into two parts, $[0, a]$ and $[a, N_0^{-1}]$ with $a$ given by:
"""

# ╔═╡ dbf0f411-344f-4d6c-b031-b5c61022cbd4
a = Arf(1 // 1024)

# ╔═╡ d8edc54e-6bfa-4640-933f-f924a9ffe38b
md"""
We then verify that $G_q^{(5)}$ is positive on $[0, a]$:
"""

# ╔═╡ a14c2a9a-4561-4420-b7aa-fe386bc84ea2
G_q_d5_minimum_0_a = ArbExtras.minimum_enclosure(
    ArbExtras.derivative_function(G_q, 5),
    Arf(0),
    a,
    verbose = true,
)

# ╔═╡ 5169dc6c-8189-460a-a7c4-37bd4360f7d3
Arblib.ispositive(G_q_d5_minimum_0_a)

# ╔═╡ 7a199b36-204a-4a9b-9a34-285b52b56741
md"""
Finally we directly enclose the minimum of $G_q$ on the interval $[a, N_0^{-1}]$ and check that it is positive:
"""

# ╔═╡ 94b731fb-019c-43e1-b0b3-519d0f801529
G_q_minimum_a_inv_N₀ =
    ArbExtras.minimum_enclosure(G_q, a, ubound(Arb(1 // N₀)), verbose = true)

# ╔═╡ cdeb73d6-600e-49b4-8dd9-2329d4b9dc9d
Arblib.ispositive(G_q_minimum_a_inv_N₀)

# ╔═╡ Cell order:
# ╟─bfb90022-f09a-46a0-b0aa-5d1ffdd7a908
# ╠═a8ef131c-85bb-11f0-28cf-7b4c670ed13b
# ╟─d2d34a83-67a8-4d71-8b06-171fefb55487
# ╠═1e814a20-d22d-4c0d-9003-6712a0a5d205
# ╟─ecdfb51c-bd88-4d2b-86fd-d6c46340bf28
# ╠═03f4726b-2d20-4a77-83ba-41824798533f
# ╟─bb882554-bd78-4ab2-acfc-e7fb702e0024
# ╠═76af3933-de21-45f0-9310-028eaf31465f
# ╟─4eb02435-b201-4e16-9ac5-fbc68c131dad
# ╠═b20f3e75-b475-4cac-b69d-b6611ac82363
# ╟─bfcf8109-a214-4460-b41a-35a714096f05
# ╠═0d216fba-d53a-47ba-886f-dc8902b0b0c2
# ╟─e03b99aa-2802-470e-adb3-dbf6a74da074
# ╠═f7c11696-e359-48c6-bf29-b799b3a26957
# ╟─b2b04733-b22f-4283-b657-013bd02864ac
# ╠═1b5615d5-ff94-41c6-8123-d11bdae8af08
# ╟─8184b4d2-7235-4cfe-8a34-a267ae170be7
# ╟─4371e67c-6c2e-461c-bbb5-afcb8ad7f142
# ╠═89947765-e753-4750-b6cc-084c0c86280f
# ╟─80973d1e-71af-454e-b0ee-a949e445ed44
# ╠═c427a061-0bf9-4e40-947e-2bd093ad686a
# ╟─40f8bc3e-9d69-4e95-800c-7f33e5187b8c
# ╠═8374aeda-b25a-4a9a-b392-47b423892b86
# ╟─844d9803-36ca-40a2-8dbf-911a52ede8dc
# ╠═218d06e7-c755-41b3-bfaa-1dfd6a3ccded
# ╟─bcbcb115-b950-4e71-ab50-665e67b32c68
# ╠═dbf0f411-344f-4d6c-b031-b5c61022cbd4
# ╟─d8edc54e-6bfa-4640-933f-f924a9ffe38b
# ╠═a14c2a9a-4561-4420-b7aa-fe386bc84ea2
# ╠═5169dc6c-8189-460a-a7c4-37bd4360f7d3
# ╟─7a199b36-204a-4a9b-9a34-285b52b56741
# ╠═94b731fb-019c-43e1-b0b3-519d0f801529
# ╠═cdeb73d6-600e-49b4-8dd9-2329d4b9dc9d
