### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 5cce8f4c-b026-11f0-0c40-9df9d697234e
begin
    using Pkg
    Pkg.activate("..", io = devnull)
    using SpectralRegularPolygons
    using Arblib

    import SpectralRegularPolygons as SRP

    setprecision(Arb, 128)
end

# ╔═╡ fd460f81-546d-4dcc-b175-477c711faecd
md"""
# Proof of Lemma C.3
"""

# ╔═╡ 4af55273-a0cf-4482-898e-6e4756c165b3
md"""
## Goal
We want to prove that for $N \geq N_0$, the $L_j$'s from Lemma C.2 have the bound $L_j \leq C_{L_j}$. Here $N_0$ and the $C_{L_j}$'s are given by
"""

# ╔═╡ ab492e70-f2bb-462a-9a06-830421e78c13
N₀ = SRP.N₀

# ╔═╡ f83262d6-3d38-418c-b771-0efab558b2a9
C_L_0 = SRP.C_L_0

# ╔═╡ b144c785-b67e-491a-9a13-92cbea1b46fa
C_L_1 = SRP.C_L_1

# ╔═╡ da042862-22af-471f-97c5-25c588b80b3d
C_L_2 = SRP.C_L_2

# ╔═╡ 901bf08f-6244-4b37-928b-0b8b57322a94
C_L_3 = SRP.C_L_3

# ╔═╡ 9157dcc1-1dbc-48d6-b774-7a1ff9e0e027
C_L_4 = SRP.C_L_4

# ╔═╡ 04dda97f-5a41-4625-b3ed-cf297b518304
C_L_6 = SRP.C_L_6

# ╔═╡ d1dea2a8-a403-41f9-9a6e-0017d8b113e4
md"""
## Proof
We have the following values for the constants coming from Lemma C.1
"""

# ╔═╡ 791fff4f-3fbd-43e8-aff0-33e3fac72d15
C_ρ_λ = Arb(SRP.C_ρ_λ)

# ╔═╡ 0a6e909e-52e5-4eb9-87c2-134745b41f7b
C_S_2 = Arb(SRP.C_S_2)

# ╔═╡ 322b4704-d5a7-4af7-b284-4f76ea9abbbf
C_S_3 = Arb(SRP.C_S_3)

# ╔═╡ 9cec456a-7b8b-4267-9207-f898cd841f36
C_F_N_0 = Arb(SRP.C_F_N_0)

# ╔═╡ 37c50b4b-72f2-4988-83c3-b64fa297a21a
C_F_N_2 = Arb(SRP.C_F_N_2)

# ╔═╡ 05c4a227-693d-4322-bc0e-61ae7048f53a
C_F_N_3 = Arb(SRP.C_F_N_3)

# ╔═╡ 6e2fd208-9320-4230-bc6c-8aeac2bc56e4
C_F_N_4 = Arb(SRP.C_F_N_4)

# ╔═╡ 9d0b358d-b083-478b-bb97-0954d878e672
C_J0_8 = Arb(SRP.C_J0_8)

# ╔═╡ b0cbad8b-e994-497a-93a0-d84ea2aedca9
md"""
In addition to the above bounds we also need an enclosure for $\lambda$:
"""

# ╔═╡ 027a7f37-2ea6-471c-9251-ce073061e6c4
λ = SRP.λ_disc()

# ╔═╡ bc51790b-b19b-476e-94d5-8e144c103634
md"""
For bounding the $L_j$'s we start by bounding the terms $T_{K,k,j}$ for $k = 1, 2, 3$. Any division by $N$ in these expressions we replace by $N_0$, since that still gives an upper bound valid for $N \geq N_0$.
"""

# ╔═╡ 71bed7dc-a8ff-4c57-b414-a1b5e927c88c
md"""
### Bounds for $T_{K,1}$
"""

# ╔═╡ bd6f63cf-76d1-474f-80f5-260f006b0266
T_K_1_0 =
    1 // 2 * C_F_N_0 * C_F_N_4 +
    1 // 2 * C_F_N_2 * (C_S_2 + C_S_3 / N₀) +
    1 // 2 * C_ρ_λ / N₀^3 * (C_S_2 + C_S_3 / N₀)

# ╔═╡ ee9031e4-4bad-4c44-81e6-6297561543e9
T_K_1_1 =
    1 // 4 * ((1 + C_F_N_0 + C_S_2 / N₀^2) * C_F_N_3) +
    1 // 2N₀ * C_F_N_2 * C_S_2 +
    1 // 4N₀ * C_S_2^2

# ╔═╡ 28b2d3fd-93d0-4fae-a31f-b71b2e6d7fec
T_K_1_2 = Arb(0)

# ╔═╡ 751cb540-1243-4f9a-a54e-6669b987944e
T_K_1_3 = 1 // 8 * (1 + C_F_N_0) * C_F_N_2

# ╔═╡ 670ead62-5742-4d68-aa7d-cc7a573d4e3f
T_K_1_4 =
    1 // 96 * C_F_N_0 * (1 + C_F_N_4 / N₀^4) +
    1 // 96N₀^2 * C_F_N_0 * (C_S_2 + C_S_3 / N₀) +
    1 // 24N₀ * (1 + C_F_N_0) * C_F_N_2

# ╔═╡ 791273d4-d3f5-4a49-898e-6e1911153108
md"""
### Bounds for $T_{K,2}$
"""

# ╔═╡ 42ff5795-514b-44dd-9710-aa95e6979b23
T_K_2_0 = 1 // 16 * C_F_N_0^2 * C_F_N_2^2

# ╔═╡ fe84ab1f-1015-4709-b28e-e77668de66b4
T_K_2_1 =
    1 // 16 * C_F_N_0^3 * C_F_N_3 +
    1 // 16N₀^4 * C_ρ_λ * C_S_2 +
    1 // 32N₀ * (1 + C_F_N_0 * (1 + C_F_N_0)) * C_F_N_2 * C_S_2


# ╔═╡ 265a5374-dda9-4368-8e62-4bcb7d2e4138
T_K_2_2 =
    1 // 32 * C_F_N_0^3 * C_F_N_2 + 1 // 64 * (1 + C_F_N_0^2) * (1 + C_F_N_0) * C_F_N_2

# ╔═╡ bf9739e9-e9a9-445f-953f-030c14e47f76
T_K_2_3 = 1 // 64N₀ * (1 + C_F_N_0^2) * (1 + C_F_N_0) * C_F_N_2

# ╔═╡ 6e750e74-219b-4af5-823d-a46e38e82b7b
T_K_2_4 =
    1 // 256 * C_F_N_0^4 +
    1 // 192 * C_F_N_0^3 * (C_F_N_0 + C_S_2 / N₀^2) +
    1 / 192N₀^4 * (1 + C_F_N_0 * (1 + C_F_N_0)) * C_F_N_2 * C_S_2 +
    1 / 192N₀^2 * C_S_2

# ╔═╡ e0a4ea5f-80d2-41c5-9090-88ac3f715256
md"""
### Bounds for $T_{K,3}$
"""

# ╔═╡ f9a81055-b26d-4190-8090-194c39b06422
T_K_3_0 = 1 // 288N₀^2 * C_F_N_0^3 * C_F_N_2^3

# ╔═╡ cc96ba5f-bb64-4531-bdd7-71c8d088efe6
T_K_3_1 = 1 // 192N₀ * C_F_N_0^4 * C_F_N_2^2

# ╔═╡ e5d45b1a-4137-4427-9c1a-583fd8de4512
T_K_3_2 = 1 // 384N₀^2 * C_F_N_0^4 * C_F_N_2^2 + 1 // 384 * C_F_N_0^5 * C_F_N_2

# ╔═╡ 20c3a605-6e62-42cd-ba98-6ef345227c58
T_K_3_3 =
    1 // 2304N₀ * (1 + C_F_N_0^3) * (1 + C_F_N_0 + C_F_N_0^2) * C_F_N_2 +
    1 // 768N₀ * C_F_N_0^5 * C_F_N_2

# ╔═╡ 497b0c84-79bd-40ed-93a7-dd1bd1f869b6
T_K_3_4 = 1 // 1536N₀^2 * C_F_N_0^5 * C_F_N_2 + 1 / 1536 * C_F_N_0^6

# ╔═╡ 5919d4c1-48d8-4267-bf08-b3db4f2668ec
T_K_3_6 = 1 / 18432N₀^2 * C_F_N_0^6

# ╔═╡ dfab2abc-32e2-4330-a1f3-57af8f9c6d52
md"""
### Bounds for $L_j$'s
We combine the above bounds with the bound for the remaning term to get bounds for the $L_j$'s.
"""

# ╔═╡ 652e1a8d-4bea-4692-9585-fd0c6b0c0fb5
L_0 =
    λ * T_K_1_0 +
    λ^2 * T_K_2_0 +
    λ^3 * T_K_3_0 +
    λ^4 * 16 // N₀^4 * C_J0_8 / factorial(8) * C_F_N_0^4 * C_F_N_4^4

# ╔═╡ aa3821e1-7e24-448d-a45c-fe38f7cc1ad8
L_1 =
    λ * T_K_1_1 +
    λ^2 * T_K_2_1 +
    λ^3 * T_K_3_1 +
    λ^4 * 32 // N₀^3 * C_J0_8 / factorial(8) * C_F_N_0^5 * C_F_N_4^3

# ╔═╡ b2e037d6-af3e-42d1-9841-d5f2515396ca
L_2 =
    λ * T_K_1_2 +
    λ^2 * T_K_2_2 +
    λ^3 * T_K_3_2 +
    λ^4 * 24 // N₀^2 * C_J0_8 / factorial(8) * C_F_N_0^6 * C_F_N_4^2

# ╔═╡ 9088108e-8d30-4f0a-b684-ed0a4ca3ddc2
L_3 =
    λ * T_K_1_3 +
    λ^2 * T_K_2_3 +
    λ^3 * T_K_3_3 +
    λ^4 * 8 // N₀ * C_J0_8 / factorial(8) * C_F_N_0^7 * C_F_N_4

# ╔═╡ 0253d16f-ed6d-48e9-af2e-88ea6071a57e
L_4 = λ * T_K_1_4 + λ^2 * T_K_2_4 + λ^3 * T_K_3_4 + λ^4 * C_J0_8 / factorial(8) * C_F_N_0^8

# ╔═╡ 49cb3d73-a689-42f3-8ec5-a4fdd7ae2a7d
L_6 = λ^3 * T_K_3_6

# ╔═╡ 5007f736-184e-45f2-86d0-089a9b2de3cd
md"""
Finally we verify that they satisfy the proposed bounds.
"""

# ╔═╡ 914a0f47-1383-4508-a09e-219db4f0093a
@assert_proof L_0 < Arb(C_L_0)

# ╔═╡ f6cbd151-c8b9-4db0-9d3c-55253b4145c8
@assert_proof L_1 < Arb(C_L_1)

# ╔═╡ bfec4b1c-2331-4758-996c-c20b994fe295
@assert_proof L_2 < Arb(C_L_2)

# ╔═╡ a92f2e6a-0509-47f0-9909-7c311863f6d9
@assert_proof L_3 < Arb(C_L_3)

# ╔═╡ 1a5f79b7-7381-434b-ae9f-573f26989e0b
@assert_proof L_4 < Arb(C_L_4)

# ╔═╡ 7f7f7b23-ac11-4d40-9600-b64f8be6b8bf
@assert_proof L_6 < Arb(C_L_6)

# ╔═╡ Cell order:
# ╟─fd460f81-546d-4dcc-b175-477c711faecd
# ╠═5cce8f4c-b026-11f0-0c40-9df9d697234e
# ╟─4af55273-a0cf-4482-898e-6e4756c165b3
# ╠═ab492e70-f2bb-462a-9a06-830421e78c13
# ╠═f83262d6-3d38-418c-b771-0efab558b2a9
# ╠═b144c785-b67e-491a-9a13-92cbea1b46fa
# ╠═da042862-22af-471f-97c5-25c588b80b3d
# ╠═901bf08f-6244-4b37-928b-0b8b57322a94
# ╠═9157dcc1-1dbc-48d6-b774-7a1ff9e0e027
# ╠═04dda97f-5a41-4625-b3ed-cf297b518304
# ╟─d1dea2a8-a403-41f9-9a6e-0017d8b113e4
# ╠═791fff4f-3fbd-43e8-aff0-33e3fac72d15
# ╠═0a6e909e-52e5-4eb9-87c2-134745b41f7b
# ╠═322b4704-d5a7-4af7-b284-4f76ea9abbbf
# ╠═9cec456a-7b8b-4267-9207-f898cd841f36
# ╠═37c50b4b-72f2-4988-83c3-b64fa297a21a
# ╠═05c4a227-693d-4322-bc0e-61ae7048f53a
# ╠═6e2fd208-9320-4230-bc6c-8aeac2bc56e4
# ╠═9d0b358d-b083-478b-bb97-0954d878e672
# ╟─b0cbad8b-e994-497a-93a0-d84ea2aedca9
# ╠═027a7f37-2ea6-471c-9251-ce073061e6c4
# ╟─bc51790b-b19b-476e-94d5-8e144c103634
# ╟─71bed7dc-a8ff-4c57-b414-a1b5e927c88c
# ╠═bd6f63cf-76d1-474f-80f5-260f006b0266
# ╠═ee9031e4-4bad-4c44-81e6-6297561543e9
# ╠═28b2d3fd-93d0-4fae-a31f-b71b2e6d7fec
# ╠═751cb540-1243-4f9a-a54e-6669b987944e
# ╠═670ead62-5742-4d68-aa7d-cc7a573d4e3f
# ╟─791273d4-d3f5-4a49-898e-6e1911153108
# ╠═42ff5795-514b-44dd-9710-aa95e6979b23
# ╠═fe84ab1f-1015-4709-b28e-e77668de66b4
# ╠═265a5374-dda9-4368-8e62-4bcb7d2e4138
# ╠═bf9739e9-e9a9-445f-953f-030c14e47f76
# ╠═6e750e74-219b-4af5-823d-a46e38e82b7b
# ╟─e0a4ea5f-80d2-41c5-9090-88ac3f715256
# ╠═f9a81055-b26d-4190-8090-194c39b06422
# ╠═cc96ba5f-bb64-4531-bdd7-71c8d088efe6
# ╠═e5d45b1a-4137-4427-9c1a-583fd8de4512
# ╠═20c3a605-6e62-42cd-ba98-6ef345227c58
# ╠═497b0c84-79bd-40ed-93a7-dd1bd1f869b6
# ╠═5919d4c1-48d8-4267-bf08-b3db4f2668ec
# ╟─dfab2abc-32e2-4330-a1f3-57af8f9c6d52
# ╠═652e1a8d-4bea-4692-9585-fd0c6b0c0fb5
# ╠═aa3821e1-7e24-448d-a45c-fe38f7cc1ad8
# ╠═b2e037d6-af3e-42d1-9841-d5f2515396ca
# ╠═9088108e-8d30-4f0a-b684-ed0a4ca3ddc2
# ╠═0253d16f-ed6d-48e9-af2e-88ea6071a57e
# ╠═49cb3d73-a689-42f3-8ec5-a4fdd7ae2a7d
# ╟─5007f736-184e-45f2-86d0-089a9b2de3cd
# ╠═914a0f47-1383-4508-a09e-219db4f0093a
# ╠═f6cbd151-c8b9-4db0-9d3c-55253b4145c8
# ╠═bfec4b1c-2331-4758-996c-c20b994fe295
# ╠═a92f2e6a-0509-47f0-9909-7c311863f6d9
# ╠═1a5f79b7-7381-434b-ae9f-573f26989e0b
# ╠═7f7f7b23-ac11-4d40-9600-b64f8be6b8bf
