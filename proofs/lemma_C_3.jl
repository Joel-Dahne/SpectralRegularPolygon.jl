### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 5cce8f4c-b026-11f0-0c40-9df9d697234e
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
N₀ = 64

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

# ╔═╡ d7d66ff9-e8a1-47a8-a0e6-430381613b09
C_S_2_tilde = Arb(SRP.C_S_2_tilde)

# ╔═╡ 42a82761-2403-4373-a59b-6e88ae832121
C_S_3_tilde = Arb(SRP.C_S_3_tilde)

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
In addition to the above bounds we also need bounds for the following values that occurr in the expressions for the $L_j$'s.
"""

# ╔═╡ 027a7f37-2ea6-471c-9251-ce073061e6c4
λ = SRP.λ_disc()

# ╔═╡ 765a89f1-0a26-4785-8cf2-9a380a824c3d
ρ = let inv_N = Arb((0, 1 // N₀))
    SRP.c_N(inv_N)^2 * SRP.λ_app(inv_N)
end

# ╔═╡ bc51790b-b19b-476e-94d5-8e144c103634
md"""
For bounding the $L_j$'s we start by bounding the terms $T_{K,k,j}$ for $k = 1, 2, 3, 4$. Any division by $N$ in these expressions we replace by $N_0$, since that still gives an upper bound valid for $N \geq N_0$.
"""

# ╔═╡ 71bed7dc-a8ff-4c57-b414-a1b5e927c88c
md"""
### Bounds for $T_{K,1}$
"""

# ╔═╡ bd6f63cf-76d1-474f-80f5-260f006b0266
T_K_1_0 =
    λ / 4 * (1 + C_F_N_2 / N₀^2) * 2C_F_N_4 +
    C_F_N_2 * λ / 4 * (C_S_2_tilde + C_S_3_tilde / N₀ + (C_S_2_tilde + C_S_3_tilde / N₀)) +
    λ / 2 * C_ρ_λ / N₀^3 * (C_S_2 + C_S_3 / N₀)

# ╔═╡ ee9031e4-4bad-4c44-81e6-6297561543e9
T_K_1_1 =
    ((2 + (C_F_N_2 + C_S_2_tilde) / N₀^2) * C_F_N_3) * λ / 4 +
    2C_F_N_2 * C_S_2 * λ / 4N₀ +
    C_S_2^2 * λ / 4N₀

# ╔═╡ 28b2d3fd-93d0-4fae-a31f-b71b2e6d7fec
T_K_1_2 = Arb(0)

# ╔═╡ 751cb540-1243-4f9a-a54e-6669b987944e
T_K_1_3 = ((2 + C_F_N_2 / N₀) * C_F_N_2) * λ / 8

# ╔═╡ 670ead62-5742-4d68-aa7d-cc7a573d4e3f
T_K_1_4 =
    λ / 4 * (1 + C_F_N_2 / N₀^2) / 24 * (1 + C_F_N_4 / N₀^4) +
    C_F_N_2 * λ / 4 * (C_S_2_tilde + C_S_3_tilde / N₀) / 24N₀^4 +
    λ / 4 * (C_S_2_tilde + C_S_3_tilde / N₀) / 24N₀^2 +
    ((2 + C_F_N_2 / N₀^2) * C_F_N_2) * λ / 24N₀

# ╔═╡ 791273d4-d3f5-4a49-898e-6e1911153108
md"""
### Bounds for $T_{K,2}$
"""

# ╔═╡ 42ff5795-514b-44dd-9710-aa95e6979b23
T_K_2_0 = λ^2 / 16 * (1 + C_F_N_2 / N₀^2) * C_F_N_2^2

# ╔═╡ fe84ab1f-1015-4709-b28e-e77668de66b4
T_K_2_1 =
    λ^2 / 16 * (1 + C_F_N_2 / N₀^2) * C_F_N_3 +
    λ^2 / 16 * C_ρ_λ * C_S_2 / N₀^4 +
    λ^2 / 32N₀ * ((1 + C_F_N_2 / N₀^2) * (2 + C_F_N_2 / N₀^2) + 1) * C_F_N_2 * C_S_2_tilde


# ╔═╡ 265a5374-dda9-4368-8e62-4bcb7d2e4138
T_K_2_2 =
    λ^2 / 32 * (1 + C_F_N_2 / N₀^2)^3 * C_F_N_2 +
    λ^2 / 64 * ((1 + C_F_N_2 / N₀^2)^2 + 1) * (2 + C_F_N_2 / N₀^2) * C_F_N_2

# ╔═╡ bf9739e9-e9a9-445f-953f-030c14e47f76
T_K_2_3 = λ^2 / 64N₀ * ((1 + C_F_N_2 / N₀^2)^2 + 1) * (2 + C_F_N_2 / N₀^2) * C_F_N_2

# ╔═╡ 6e750e74-219b-4af5-823d-a46e38e82b7b
T_K_2_4 =
    λ^2 / 256 * (1 + C_F_N_2 / N₀^2)^4 +
    λ^2 / 192 * (1 + C_F_N_2 / N₀^2)^3 * (1 + (C_F_N_2 + C_S_2_tilde) / N₀^2) +
    λ^2 / 192N₀^4 *
    ((1 + C_F_N_2 / N₀^2) * (2 + C_F_N_2 / N₀^2) + 1) *
    C_F_N_2 *
    C_S_2_tilde +
    λ^2 / 192N₀^2 * C_S_2_tilde

# ╔═╡ e0a4ea5f-80d2-41c5-9090-88ac3f715256
md"""
### Bounds for $T_{K,3}$
"""

# ╔═╡ f9a81055-b26d-4190-8090-194c39b06422
T_K_3_0 = λ^3 / 2304 * (1 + C_F_N_2 / N₀^2)^3 * 8C_F_N_2^3 / N₀^2

# ╔═╡ cc96ba5f-bb64-4531-bdd7-71c8d088efe6
T_K_3_1 = 3λ^3 / 2304 * (1 + C_F_N_2 / N₀^2)^4 * 4C_F_N_2^2 / N₀

# ╔═╡ e5d45b1a-4137-4427-9c1a-583fd8de4512
T_K_3_2 =
    λ^3 / 2304 * (1 + C_F_N_2 / N₀^2)^3 * 12C_F_N_2^2 / 2 * (1 + C_F_N_2 / N₀^2) / N₀^2 +
    3λ^3 / 2304 * (1 + C_F_N_2 / N₀^2)^5 * 2C_F_N_2

# ╔═╡ 20c3a605-6e62-42cd-ba98-6ef345227c58
T_K_3_3 =
    λ^3 / 2304 *
    ((1 + C_F_N_2 / N₀^2)^3 + 1) *
    ((1 + C_F_N_2 / N₀^2)^2 + (2 + C_F_N_2 / N₀^2)) *
    C_F_N_2 + 3λ^3 / 2304 * (1 + C_F_N_2 / N₀^2)^5 * C_F_N_2 / N₀

# ╔═╡ 497b0c84-79bd-40ed-93a7-dd1bd1f869b6
T_K_3_4 =
    λ^3 / 2304 * (1 + C_F_N_2 / N₀^2)^3 * 6C_F_N_2 / 4 * (1 + C_F_N_2 / N₀^2)^2 / N₀^2 +
    3λ^3 / 2304 * (1 + C_F_N_2 / N₀^2)^6 / 2

# ╔═╡ 5919d4c1-48d8-4267-bf08-b3db4f2668ec
T_K_3_6 = λ^3 / 2304 * (1 + C_F_N_2 / N₀^2)^6 / 8N₀^2

# ╔═╡ b276a018-77ff-484d-83e4-c13391d49562
md"""
### Bounds for $T_{K,4}$
"""

# ╔═╡ 1ad65ea8-e6d2-4e89-b01b-60dc57ddea82
T_K_4_0 = C_J0_8 / factorial(8) * ρ^4 * (1 + C_F_N_2 / N₀^2)^4 * 16C_F_N_4^4 / N₀^4

# ╔═╡ 276b29ec-0993-45ad-a2e6-8cf5138431ba
T_K_4_1 = 4C_J0_8 / factorial(8) * ρ^4 * (1 + C_F_N_2 / N₀^2)^5 * 8C_F_N_4^3 / N₀^3

# ╔═╡ ed385d3d-2b44-465f-80ab-b246aef14dc4
T_K_4_2 = 6C_J0_8 / factorial(8) * ρ^4 * (1 + C_F_N_2 / N₀^2)^6 * 4C_F_N_4^2 / N₀^2

# ╔═╡ f0887698-38b9-4d16-ba8f-d1432862d9b9
T_K_4_3 = 4C_J0_8 / factorial(8) * ρ^4 * (1 + C_F_N_2 / N₀^2)^7 * 2C_F_N_4 / N₀

# ╔═╡ 4ae121b4-85e1-43cd-bc08-b5caf03dbd30
T_K_4_4 = C_J0_8 / factorial(8) * ρ^4 * (1 + C_F_N_2 / N₀^2)^8

# ╔═╡ dfab2abc-32e2-4330-a1f3-57af8f9c6d52
md"""
### Bounds for $L_j$'s
We combine the above bounds to get bounds for the $L_j$'s.
"""

# ╔═╡ 652e1a8d-4bea-4692-9585-fd0c6b0c0fb5
L_0 = T_K_1_0 + T_K_2_0 + T_K_3_0 + T_K_4_0

# ╔═╡ aa3821e1-7e24-448d-a45c-fe38f7cc1ad8
L_1 = T_K_1_1 + T_K_2_1 + T_K_3_1 + T_K_4_1

# ╔═╡ b2e037d6-af3e-42d1-9841-d5f2515396ca
L_2 = T_K_1_2 + T_K_2_2 + T_K_3_2 + T_K_4_2

# ╔═╡ 9088108e-8d30-4f0a-b684-ed0a4ca3ddc2
L_3 = T_K_1_3 + T_K_2_3 + T_K_3_3 + T_K_4_3

# ╔═╡ 0253d16f-ed6d-48e9-af2e-88ea6071a57e
L_4 = T_K_1_4 + T_K_2_4 + T_K_3_4 + T_K_4_4

# ╔═╡ 49cb3d73-a689-42f3-8ec5-a4fdd7ae2a7d
L_6 = T_K_3_6

# ╔═╡ 5007f736-184e-45f2-86d0-089a9b2de3cd
md"""
Finally we verify that they satisfy the proposed bounds.
"""

# ╔═╡ 914a0f47-1383-4508-a09e-219db4f0093a
L_0 < Arb(C_L_0)

# ╔═╡ f6cbd151-c8b9-4db0-9d3c-55253b4145c8
L_1 < Arb(C_L_1)

# ╔═╡ bfec4b1c-2331-4758-996c-c20b994fe295
L_2 < Arb(C_L_2)

# ╔═╡ a92f2e6a-0509-47f0-9909-7c311863f6d9
L_3 < Arb(C_L_3)

# ╔═╡ 1a5f79b7-7381-434b-ae9f-573f26989e0b
L_4 < Arb(C_L_4)

# ╔═╡ 7f7f7b23-ac11-4d40-9600-b64f8be6b8bf
L_6 < Arb(C_L_6)

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
# ╠═d7d66ff9-e8a1-47a8-a0e6-430381613b09
# ╠═42a82761-2403-4373-a59b-6e88ae832121
# ╠═9cec456a-7b8b-4267-9207-f898cd841f36
# ╠═37c50b4b-72f2-4988-83c3-b64fa297a21a
# ╠═05c4a227-693d-4322-bc0e-61ae7048f53a
# ╠═6e2fd208-9320-4230-bc6c-8aeac2bc56e4
# ╠═9d0b358d-b083-478b-bb97-0954d878e672
# ╟─b0cbad8b-e994-497a-93a0-d84ea2aedca9
# ╠═027a7f37-2ea6-471c-9251-ce073061e6c4
# ╠═765a89f1-0a26-4785-8cf2-9a380a824c3d
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
# ╟─b276a018-77ff-484d-83e4-c13391d49562
# ╠═1ad65ea8-e6d2-4e89-b01b-60dc57ddea82
# ╠═276b29ec-0993-45ad-a2e6-8cf5138431ba
# ╠═ed385d3d-2b44-465f-80ab-b246aef14dc4
# ╠═f0887698-38b9-4d16-ba8f-d1432862d9b9
# ╠═4ae121b4-85e1-43cd-bc08-b5caf03dbd30
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
