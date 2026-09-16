### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 5ccffc41-fd3b-439c-86bb-ec92e86fa30a
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

# ╔═╡ 9047fbda-85bb-11f0-3c79-517d4f9a6381
md"""
# Proof of Lemma 2.21
"""

# ╔═╡ a8e5d438-7f29-43c5-8576-eccff74668d9
md"""
## Goal
We want to prove that for $N \geq N_0$ we have

$$\frac{\lambda_{\text{app}}}{1 - \hat{\varepsilon}(N)} < \lambda_2(\mathbb{D}_{R_{\text{outer}}}),$$

where $\lambda_2(\mathbb{D}_{R_{\text{outer}}})$ denotes the second eigenvalue of the $\mathbb{D}_{R_{\text{outer}}}$ and $N_0$ and $R_{\text{outer}}$ are given by:
"""

# ╔═╡ 6a32743d-ec5d-44be-9b84-20c55481446f
N₀ = SRP.N₀

# ╔═╡ 8c1beb93-8fba-4269-9001-204575c784f1
R_outer = SRP.R_outer

# ╔═╡ 79db37ba-4245-47ee-95fb-3d9aa2a48e40
md"""
## Proof
### Step 1 - Compute $\lambda_2(\mathbb{D}_{R_{\text{outer}}})$
As a first step we compute an enclosure of $\lambda_2(\mathbb{D}_{R_{\text{outer}}})$. By scaling we get that

$$\lambda_2(\mathbb{D}_{R_{\text{outer}}}) = \frac{\lambda_2(\mathbb{D})}{R_{\text{outer}}^2}.$$

To compute $\lambda_2(\mathbb{D})$ we use that it is given by

$$\lambda_2(\mathbb{D}) = j_{1,1}^2,$$

where $j_{1,1}$ denotes the first positive zero of $J_1$.
"""

# ╔═╡ f014ea10-6ce2-409c-aceb-8191e0debe4f
md"""
To compute $j_{1,1}$ we first verify that $J_1$ is strictly increasing on the interval $[0, 1]$, this ensures that the only zero on that interval is the one at zero.
"""

# ╔═╡ fe4a5de7-8ba3-42f3-ae35-19de44b4ded3
Arblib.ispositive(ArbExtras.derivative_function(besselj1)(Arb((0, 1))))

# ╔═╡ 725399ac-9ed3-4e92-a38c-725668df4328
md"""
This proves that $j_{1,1} > 1$. Next we isolate all the roots on the interval $[1, 5]$:
"""

# ╔═╡ ccb58a39-cdf1-45d3-9dc0-8eabe65021b8
roots, flags = ArbExtras.isolate_roots(besselj1, Arf(1), Arf(5), depth = 20)

# ╔═╡ fc8afe74-e745-487b-9303-98b9341e665c
md"""
We verify that it only found one root and that it was proved to be unique:
"""

# ╔═╡ 6b68bc60-829d-4a32-8bc7-33644b6b1edc
length(flags) == 1 && flags[1]

# ╔═╡ 07199fe9-c599-4a9e-8066-538d9c56b95c
md"""
We refine the enclosure of the root, giving an enclosure for $j_{1,1}$:
"""

# ╔═╡ 0aaaf5f1-8728-4648-b09d-250376784390
j_1_1 = ArbExtras.refine_root(besselj1, Arb(roots[1]))

# ╔═╡ 2943ea8e-626d-43c1-ad8b-3c704b7328b0
md"""
We square it to get an enclosure of the second eigenvalue of $\mathbb{D}$:
"""

# ╔═╡ 9ff56c39-723d-4206-8895-b11bb19ab4f7
λ₂ = j_1_1^2

# ╔═╡ 4550ca4a-557a-4508-991b-d9150d6bf4f2
md"""
Finally we scale the result to get an enclosure of $\lambda_2(\mathbb{D}_{R_{\text{outer}}})$:
"""

# ╔═╡ dc6c68cb-2624-449d-b396-de0c01a29275
λ₂_R_outer = λ₂ / Arb(R_outer)^2

# ╔═╡ 150348d3-840b-4dd0-abcf-d9e4459610b0
md"""
We print a version with fewer digits for inclusion in the paper.
"""

# ╔═╡ 1863ce44-3b91-46e5-be48-00696a001688
string(λ₂_R_outer, digits = 5)

# ╔═╡ f225159f-9f55-46c0-a25c-f079eb69fe51
md"""
### Step 2 - Compute $\lambda_{\text{app}}$ and $\hat{\varepsilon}(N)$
The functions for computing $\lambda_{\text{app}}$ and $\hat{\varepsilon}(N)$ are implemented in the package. They take as input an enclosure of $N^{-1}$, so we start by computing an enclosure of this that is valid for all $N \geq N_0$:
"""

# ╔═╡ ea12574f-90c9-43d1-93f4-847f5b52351a
inv_N = Arb((0, 1 // N₀))

# ╔═╡ d709f4f2-8ca5-43b7-9087-45ec56b041bc
md"""
The values can now be enclosed:
"""

# ╔═╡ 1c759096-6aa3-4e81-bd4b-4ab753bb73f8
λ_app = SRP.λ_app(inv_N)

# ╔═╡ 0566b7a3-f43d-4146-bc69-ac541e18ea2f
ε_hat = SRP.epsilon_hat(inv_N)

# ╔═╡ 27563f0b-4c28-408e-9db2-d3f756762d6b
md"""
This gives us:
"""

# ╔═╡ bd7280e8-d4d4-4230-8909-912c3f0ae084
λ_app / (1 - ε_hat)

# ╔═╡ a10e7f9e-33c5-4f8d-b28f-587dfd9063d7
md"""
### Step 3 - Verify inequality
Finally, we just verify the inequality:
"""

# ╔═╡ 840d2b7b-8e6d-44ce-a8d5-dfa45095be07
@assert_proof λ_app / (1 - ε_hat) < λ₂_R_outer

# ╔═╡ Cell order:
# ╟─9047fbda-85bb-11f0-3c79-517d4f9a6381
# ╠═5ccffc41-fd3b-439c-86bb-ec92e86fa30a
# ╟─a8e5d438-7f29-43c5-8576-eccff74668d9
# ╠═6a32743d-ec5d-44be-9b84-20c55481446f
# ╠═8c1beb93-8fba-4269-9001-204575c784f1
# ╟─79db37ba-4245-47ee-95fb-3d9aa2a48e40
# ╟─f014ea10-6ce2-409c-aceb-8191e0debe4f
# ╠═fe4a5de7-8ba3-42f3-ae35-19de44b4ded3
# ╟─725399ac-9ed3-4e92-a38c-725668df4328
# ╠═ccb58a39-cdf1-45d3-9dc0-8eabe65021b8
# ╟─fc8afe74-e745-487b-9303-98b9341e665c
# ╠═6b68bc60-829d-4a32-8bc7-33644b6b1edc
# ╟─07199fe9-c599-4a9e-8066-538d9c56b95c
# ╠═0aaaf5f1-8728-4648-b09d-250376784390
# ╟─2943ea8e-626d-43c1-ad8b-3c704b7328b0
# ╠═9ff56c39-723d-4206-8895-b11bb19ab4f7
# ╟─4550ca4a-557a-4508-991b-d9150d6bf4f2
# ╠═dc6c68cb-2624-449d-b396-de0c01a29275
# ╟─150348d3-840b-4dd0-abcf-d9e4459610b0
# ╠═1863ce44-3b91-46e5-be48-00696a001688
# ╟─f225159f-9f55-46c0-a25c-f079eb69fe51
# ╠═ea12574f-90c9-43d1-93f4-847f5b52351a
# ╟─d709f4f2-8ca5-43b7-9087-45ec56b041bc
# ╠═1c759096-6aa3-4e81-bd4b-4ab753bb73f8
# ╠═0566b7a3-f43d-4146-bc69-ac541e18ea2f
# ╟─27563f0b-4c28-408e-9db2-d3f756762d6b
# ╠═bd7280e8-d4d4-4230-8909-912c3f0ae084
# ╟─a10e7f9e-33c5-4f8d-b28f-587dfd9063d7
# ╠═840d2b7b-8e6d-44ce-a8d5-dfa45095be07
