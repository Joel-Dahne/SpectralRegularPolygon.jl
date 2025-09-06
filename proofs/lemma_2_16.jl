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
    using OhMyThreads
    using PlutoUI
    using SpecialFunctions

    import SpectralRegularPolygons as SRP
    import ProgressLogging: @withprogress, @logprogress

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
and then verify that it satisfies the requirement. Note that it suffices to verify that $f$ maps $\mathbb{D}_R$ to a set that contains $\mathbb{D}_{R_{\text{inner}}}$. 

**TODO:** Explain better what we do here
"""

# ╔═╡ a71afcc0-01e6-4b7f-a9b0-86ccf4e9c9b6
res1 = map(N₀:128) do N
    let inv_N = Acb(1 // N)
        SRP._c_N(real(inv_N)) * R * real(SRP.hypgeom2f1(2inv_N, inv_N, 1 + inv_N, Acb(-R^N)))
    end
end

# ╔═╡ c0433aef-5b65-46ee-83f8-bfbccc9cd0df
all(Arb(R_inner) .< res1)

# ╔═╡ 4c40fe01-bbb6-442d-af53-9413cdb1a20f
res2 = let inv_N = Acb(Arb((0, 1 // 128)))
    SRP.c_N(N₀) * R * real(SRP.hypgeom2f1(2inv_N, inv_N, 1 + inv_N, Acb(Arb((-R^128, 0)))))
end

# ╔═╡ bcbc717a-f1e2-4acf-93b9-392abc3dfc4f
Arb(R_inner) < res2

# ╔═╡ b973d37b-9a60-4104-9274-b35a05a6ab35
md"""
### Step 2: Bound functions on $\mathbb{D}_R$
From step 1 we get that to bound

$$\left|\frac{V_l(f^{-1}(x)^N)}{f'(f^{-1}(x))}\right| \leq C_{V_l},$$

for $x \in \mathbb{D}_{R_{\text{inner}}}$ it suffices to bound

$$\left|\frac{V_l(z^N)}{f'(z)}\right| \leq C_{V_l}$$

for $z \in \mathbb{D}_R.$ Note that both numerator and denominator are non-zero, we can hence bound them individually.
"""

# ╔═╡ eeed099f-6ab6-47a4-802c-3c15c1d77acd
md"""
We have

$$f'(z) = \frac{c_N}{(1 - z^N)^{2/N}}.$$

For $N \geq N_0$ and $z \in \mathbb{D}_R$ we can thus directly compute an enclosure as:
"""

# ╔═╡ 5171eccd-39cc-4afc-9d8b-3d5c5aa17968
f_derivative = let
    inv_N = Arb((0, 1 // N₀))
    z_pow_N = add_error(Acb(0), R^N₀)
    SRP.c_N(N₀) / (1 - z_pow_N)^2inv_N
end

# ╔═╡ 0ae4fcea-8a49-456d-9f93-0e7d820e023c
md"""
For the $V$'s the enclosure can also be directly computer:
"""

# ╔═╡ 4b747119-b501-4dc8-8a11-10cd7c6a58c7
V_1, V_2, V_3, V_4 = let
    z_pow_N = add_error(Acb(0), R^N₀)
    SRP.V(1, z_pow_N), SRP.V(2, z_pow_N), SRP.V(3, z_pow_N), SRP.V(4, z_pow_N)
end

# ╔═╡ 6fc302f2-f256-4ffd-ad6d-ca28edb7172a
md"""
Finally, we can compute the quotients and verify that the bounds hold:
"""

# ╔═╡ 306843ed-caa8-46a4-b502-b9eb29f50f6a
V_1_div_f_derivative = abs(V_1 / f_derivative)

# ╔═╡ f362faec-92e8-400c-90c3-84f241e4ef42
V_2_div_f_derivative = abs(V_2 / f_derivative)

# ╔═╡ 2fe32659-bb87-42ea-aafb-fdfa19cbdbdc
V_3_div_f_derivative = abs(V_3 / f_derivative)

# ╔═╡ ee0ab4d2-2244-417e-a25a-b37913349aee
V_4_div_f_derivative = abs(V_4 / f_derivative)

# ╔═╡ cd20b05d-021c-4e31-903e-e2413dda41eb
V_1_div_f_derivative <= Arb(C_V_1)

# ╔═╡ fd4cf253-9734-497d-ba47-03f5a8aa2d45
V_2_div_f_derivative <= Arb(C_V_2)

# ╔═╡ 3b72ba1e-fd40-4317-a3d4-534bb5cd4eff
V_3_div_f_derivative <= Arb(C_V_3)

# ╔═╡ 23b5a202-6d64-4288-9b7b-870939c2370b
V_4_div_f_derivative <= Arb(C_V_4)

# ╔═╡ Cell order:
# ╟─3f5fec32-8fd2-4cf5-9c74-e380a5853da1
# ╠═580af86c-85b6-11f0-05e7-e9613df5232a
# ╠═0c0d1542-c675-4aaf-8bdf-07221c25e8eb
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
# ╠═a71afcc0-01e6-4b7f-a9b0-86ccf4e9c9b6
# ╠═c0433aef-5b65-46ee-83f8-bfbccc9cd0df
# ╠═4c40fe01-bbb6-442d-af53-9413cdb1a20f
# ╠═bcbc717a-f1e2-4acf-93b9-392abc3dfc4f
# ╟─b973d37b-9a60-4104-9274-b35a05a6ab35
# ╟─eeed099f-6ab6-47a4-802c-3c15c1d77acd
# ╠═5171eccd-39cc-4afc-9d8b-3d5c5aa17968
# ╟─0ae4fcea-8a49-456d-9f93-0e7d820e023c
# ╠═4b747119-b501-4dc8-8a11-10cd7c6a58c7
# ╟─6fc302f2-f256-4ffd-ad6d-ca28edb7172a
# ╠═306843ed-caa8-46a4-b502-b9eb29f50f6a
# ╠═f362faec-92e8-400c-90c3-84f241e4ef42
# ╠═2fe32659-bb87-42ea-aafb-fdfa19cbdbdc
# ╠═ee0ab4d2-2244-417e-a25a-b37913349aee
# ╠═cd20b05d-021c-4e31-903e-e2413dda41eb
# ╠═fd4cf253-9734-497d-ba47-03f5a8aa2d45
# ╠═3b72ba1e-fd40-4317-a3d4-534bb5cd4eff
# ╠═23b5a202-6d64-4288-9b7b-870939c2370b
