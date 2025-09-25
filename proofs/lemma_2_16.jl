### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 580af86c-85b6-11f0-05e7-e9613df5232a
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

# ╔═╡ 3f5fec32-8fd2-4cf5-9c74-e380a5853da1
md"""
# Proof of Lemma 2.16
This notebook contains the computer-assisted part of the proof of Lemma 2.16.
"""

# ╔═╡ 0c0d1542-c675-4aaf-8bdf-07221c25e8eb
md"""
## Goal
We want to prove that for $N \geq N_0$ and all $x \in \mathbb{D}_{R_{\text{inner}}}$ we have

$$\left|\frac{V_l(f^{-1}(x)^N)}{f'(f^{-1}(x))}\right| \leq C_{V_l},$$

for $l = 1,\ 2,\ 3,\ 4$ and where $C_{V_l}$ are given by:
"""

# ╔═╡ aff645f1-6bc1-4b18-878e-9bdfce4ee2ae
C_V_1 = SRP.C_V_1

# ╔═╡ c6db017d-4648-4cf1-a242-38305afca079
C_V_2 = SRP.C_V_2

# ╔═╡ b68730ae-0b72-4ddb-86f1-a5b56d55cb60
C_V_3 = SRP.C_V_3

# ╔═╡ f2ffa3dc-4cfa-431a-8902-5c01d9fd6b42
C_V_4 = SRP.C_V_4

# ╔═╡ 9c53b9ef-4870-43b8-b955-0b3cc68d542a
md"""
We also have that $N_0$ and $R_{\text{inner}}$ are given by:
"""

# ╔═╡ f4c0e02d-2ac3-4acd-9d66-addfd551fca8
N₀ = SRP.N₀

# ╔═╡ d59a52c3-4281-4964-9567-fbaad28148a4
R_inner = SRP.R_inner

# ╔═╡ 92ae3ec0-8330-4f72-b093-6cb4184c67e0
md"""
## Proof
### Step 1 - bound $f^{-1}(x)$
As a first step, we want to find an $R$ such that for all $x \in \mathbb{D}_{R_{\text{inner}}}$ we have $f^{-1}(x) \in \mathbb{D}_R$. For this we make the following guess for $R$:
"""

# ╔═╡ 7724b725-ad79-4f07-b128-92a778ffaccf
R = Arb("0.951")

# ╔═╡ e370ed08-5c81-4725-b2b1-d7ff3671859a
md"""
and then verify that it satisfies the requirement. Note that since $f$ is continuous it suffices to verify that $f$ maps $\mathbb{D}_R$ to a set that contains $\mathbb{D}_{R_{\text{inner}}}$. We therefore need to show that for $|z| = R$ we have

$$|f(z)| > R_{\text{inner}}.$$

As discussed in the paper we have the bound

$$|f(z)| \geq c_{N_0}R\left|{}_2F_1 \left(\frac{2}{N}, \frac{1}{N}, 1 + \frac{1}{N}; -R^N\right)\right|.$$

Enclosing the ${}_2 F_1$ function we get
"""

# ╔═╡ aab12f29-6a10-4a4b-aed2-b38586a11030
hypgeom2f1_lower_bound = let inv_N = Acb(Arb((0, 1 // N₀)))
    abs(SRP.hypgeom2f1(2inv_N, inv_N, 1 + inv_N, Acb(Arb((-R^N₀, 0)))))
end

# ╔═╡ 0e8957bf-ba22-452e-8298-431e755a0248
md"""
Giving the following lower bound for $f(z)$
"""

# ╔═╡ c8b728d8-bb50-485c-a07c-b94ee856134b
f_lower_bound = SRP.c_N(Arb(1 // N₀)) * R * hypgeom2f1_lower_bound

# ╔═╡ 3c7e96b8-4041-4449-84a5-926464afda10
md"""
We can now veryify that indeed $|f(z)| > R_{\text{inner}}$
"""

# ╔═╡ c0433aef-5b65-46ee-83f8-bfbccc9cd0df
f_lower_bound > Arb(R_inner)

# ╔═╡ b973d37b-9a60-4104-9274-b35a05a6ab35
md"""
### Step 2: Bound functions on $\mathbb{D}_R$
From step 1 we get that to bound

$$\left|\frac{V_l(f^{-1}(x)^N)}{f'(f^{-1}(x))}\right| \leq C_{V_l},$$

for $x \in \mathbb{D}_{R_{\text{inner}}}$ it suffices to bound

$$\left|\frac{V_l(z^N)}{f'(z)}\right| \leq C_{V_l}$$

for $z \in \mathbb{D}_R$.
"""

# ╔═╡ eeed099f-6ab6-47a4-802c-3c15c1d77acd
md"""
We have

$$f'(z) = \frac{c_N}{(1 - z^N)^{2/N}}.$$

The monotonicity of $c_N$ then gives us the lower bound

$$|f'(z)| \geq c_{N_0}.$$
"""

# ╔═╡ d27a57f5-a1a2-46c4-be8a-f8f873171474
md"""
For $V_l(z^N)$ we have the bound

$$|V_l(z^N)| \leq D_{l,R^{N_0}}R^{N_0}$$

with which we can compute the bounds
"""

# ╔═╡ 306843ed-caa8-46a4-b502-b9eb29f50f6a
V_1_div_f_derivative_bound = SRP.V_1_div_z_bound(R^N₀) * R^N₀ / SRP.c_N(Arb(1 // N₀))

# ╔═╡ 3e005f47-3652-4e42-ae1e-aa45ca1501d3
V_2_div_f_derivative_bound = SRP.V_2_div_z_bound(R^N₀) * R^N₀ / SRP.c_N(Arb(1 // N₀))

# ╔═╡ 7cf179cf-1bd1-476d-afb3-b0cd1878ce57
V_3_div_f_derivative_bound = SRP.V_3_div_z_bound(R^N₀) * R^N₀ / SRP.c_N(Arb(1 // N₀))

# ╔═╡ cfafb8aa-1b68-4eee-8879-279848a2a816
V_4_div_f_derivative_bound = SRP.V_4_div_z_bound(R^N₀) * R^N₀ / SRP.c_N(Arb(1 // N₀))

# ╔═╡ 937a9361-dbc9-4a37-965f-c67c1044e173
md"""
We can now verify that the inequalities hold.
"""

# ╔═╡ 296b6791-55d5-4e3a-a588-07e532a50233
V_1_div_f_derivative_bound <= Arb(C_V_1)

# ╔═╡ 534cae05-8e68-4e71-b9c0-da59d8d6a5a0
V_2_div_f_derivative_bound <= Arb(C_V_2)

# ╔═╡ d0f31778-4df0-421d-af93-2cb282d56428
V_3_div_f_derivative_bound <= Arb(C_V_3)

# ╔═╡ e2294887-d1f1-4685-9db6-a7a4d88575bb
V_4_div_f_derivative_bound <= Arb(C_V_4)

# ╔═╡ fc28d9db-e652-4cd4-9052-ad57e85d1735
md"""
We print less precise bounds for inclusing in the paper.
"""

# ╔═╡ 960ac5af-fb5c-4ad7-9c26-bef9aeb549b6
string(V_1_div_f_derivative_bound, digits = 5)

# ╔═╡ 20a841ac-5687-43bc-8c9a-4d5435946733
string(V_2_div_f_derivative_bound, digits = 5)

# ╔═╡ f99dfdd9-188a-46e9-9489-42b5e4f5329e
string(V_3_div_f_derivative_bound, digits = 5)

# ╔═╡ 766ded14-8337-432d-8080-91a44b90807a
string(V_4_div_f_derivative_bound, digits = 5)

# ╔═╡ Cell order:
# ╟─3f5fec32-8fd2-4cf5-9c74-e380a5853da1
# ╠═580af86c-85b6-11f0-05e7-e9613df5232a
# ╟─0c0d1542-c675-4aaf-8bdf-07221c25e8eb
# ╠═aff645f1-6bc1-4b18-878e-9bdfce4ee2ae
# ╠═c6db017d-4648-4cf1-a242-38305afca079
# ╠═b68730ae-0b72-4ddb-86f1-a5b56d55cb60
# ╠═f2ffa3dc-4cfa-431a-8902-5c01d9fd6b42
# ╟─9c53b9ef-4870-43b8-b955-0b3cc68d542a
# ╠═f4c0e02d-2ac3-4acd-9d66-addfd551fca8
# ╠═d59a52c3-4281-4964-9567-fbaad28148a4
# ╟─92ae3ec0-8330-4f72-b093-6cb4184c67e0
# ╠═7724b725-ad79-4f07-b128-92a778ffaccf
# ╟─e370ed08-5c81-4725-b2b1-d7ff3671859a
# ╠═aab12f29-6a10-4a4b-aed2-b38586a11030
# ╟─0e8957bf-ba22-452e-8298-431e755a0248
# ╠═c8b728d8-bb50-485c-a07c-b94ee856134b
# ╟─3c7e96b8-4041-4449-84a5-926464afda10
# ╠═c0433aef-5b65-46ee-83f8-bfbccc9cd0df
# ╟─b973d37b-9a60-4104-9274-b35a05a6ab35
# ╟─eeed099f-6ab6-47a4-802c-3c15c1d77acd
# ╟─d27a57f5-a1a2-46c4-be8a-f8f873171474
# ╠═306843ed-caa8-46a4-b502-b9eb29f50f6a
# ╠═3e005f47-3652-4e42-ae1e-aa45ca1501d3
# ╠═7cf179cf-1bd1-476d-afb3-b0cd1878ce57
# ╠═cfafb8aa-1b68-4eee-8879-279848a2a816
# ╟─937a9361-dbc9-4a37-965f-c67c1044e173
# ╠═296b6791-55d5-4e3a-a588-07e532a50233
# ╠═534cae05-8e68-4e71-b9c0-da59d8d6a5a0
# ╠═d0f31778-4df0-421d-af93-2cb282d56428
# ╠═e2294887-d1f1-4685-9db6-a7a4d88575bb
# ╟─fc28d9db-e652-4cd4-9052-ad57e85d1735
# ╠═960ac5af-fb5c-4ad7-9c26-bef9aeb549b6
# ╠═20a841ac-5687-43bc-8c9a-4d5435946733
# ╠═f99dfdd9-188a-46e9-9489-42b5e4f5329e
# ╠═766ded14-8337-432d-8080-91a44b90807a
