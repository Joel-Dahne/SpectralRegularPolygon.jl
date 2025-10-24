### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 9658d218-adae-11f0-1a53-675306a5a87f
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

# ╔═╡ 0823e7e3-1b4d-409e-9668-150149d8aa1e
md"""
# Proof of Lemma C.1
"""

# ╔═╡ 8c455525-a5ec-43d1-9afb-b31c8668a15e
md"""
## Goal
We want to establish a number of bounds related to the expressions in Lemma C.2. In general these bound should hold for $N \geq N_0$ and $|t| \leq 1$, with $N_0$ given by
"""

# ╔═╡ 9362bb77-b3d5-463f-9962-229f832b3497
N₀ = SRP.N₀

# ╔═╡ f4ceb785-492f-4294-983e-c29360febd5f
md"""
We split the computation of the bounds into 4 parts, each treating related bounds.
"""

# ╔═╡ b037f618-1929-4c93-a452-5606cc59bb41
md"""
## Bounding $\rho^{1/2}/\lambda^{1/2}$
In this section we want to bound

$$\left|\frac{\rho^{1/2}}{\lambda^{1/2}} - 1 \right|$$

More precisely we want to verify that

$$\left|\frac{\rho^{1/2}}{\lambda^{1/2}} - 1 \right| \leq \frac{C_{\rho,\lambda}}{N^5}.$$

with
"""

# ╔═╡ abca81b1-bfdc-4b56-bae7-3d279ec14a29
C_ρ_λ = SRP.C_ρ_λ

# ╔═╡ 4711295e-0858-4833-b8e1-12384856ecad
md"""
Recall that

$$\rho = c_N^2 \lambda_{\text{app}}$$

and by construction $\lambda_{\text{app}}$ is taken such that for the Taylor expansion at $N = \infty$ all terms up to $N^{-4}$ in 

$$\frac{\rho^{1/2}}{\lambda^{1/2}} - 1$$

vanish. It therefore suffices to compute a bound for the fifth term in the Taylor expansion valid for $N \geq N_0$. This is equivalent to enclosing the fifth derivative in $N^{-1}$ with $N^{-1} \in [0, 1 / N_0]$ and dividing it by $5!$.
"""

# ╔═╡ 626e33ce-ece6-4066-af07-98bcb6c053c8
# Function for computing fifth derivative divided by 5!
sqrt_ρ_div_λ_d5_function = ArbExtras.derivative_function(5) do inv_N
    SRP.c_N(inv_N) * sqrt(SRP.λ_app_div_λ(inv_N)) / factorial(5)
end

# ╔═╡ a57c2ae1-85f3-4940-8d9b-2f4680fd79ef
md"""
For this we get the enclosure
"""

# ╔═╡ 7c4055e6-60b9-4085-b9c8-7fe4da7c75d7
sqrt_ρ_div_λ_5 =
    ArbExtras.enclosure_series(sqrt_ρ_div_λ_d5_function, Arb((0, 1 // N₀)), degree = 4)

# ╔═╡ a317e983-b641-4ca5-99a3-2b58fd0625c5
md"""
Which we can verify satisfies the required inequality:
"""

# ╔═╡ c3c77e89-9059-45d6-a1e4-ab5f4a06d584
abs(sqrt_ρ_div_λ_5) <= Arb(C_ρ_λ)

# ╔═╡ c3264d17-3129-4cb7-a4d9-5e7acca2c2f8
md"""
## Bounding $S_2$ and $S_3$
In this section we want to bound $S_2(t)$ and $S_3(t)$ for $|t| \leq 1$. In particular we want to show that

$$|S_2(t)| \leq C_{S_2},\quad |S_3(t)| \leq C_{S_3}$$

with
"""

# ╔═╡ e6deace5-c025-482e-bc74-b8f66da238dd
C_S_2 = SRP.C_S_2

# ╔═╡ 4b113a97-293f-4644-87fd-f536014caa71
C_S_3 = SRP.C_S_3

# ╔═╡ 8a899965-2cce-4ce2-aca1-a9953a5ad9a8
md"""
**TODO:** Note that $\frac{\rho^{1/2}}{\lambda^{1/2}}$ is bounded by $1$, so we can probably remove this in the end.

Moreover we want to show that

$$\left|\frac{\rho^{1/2}}{\lambda^{1/2}}S_2(t)\right| \leq \tilde{C}_{S_2},\quad \left|\frac{\rho^{1/2}}{\lambda^{1/2}}S_3(t)\right| \leq \tilde{C}_{S_3}$$
with
"""

# ╔═╡ 1f71ede2-30f8-4eb5-9faf-00137e2d64e8
C_S_2_tilde = SRP.C_S_2_tilde

# ╔═╡ c97d4288-2112-48ed-889b-b00aaf4d925d
C_S_3_tilde = SRP.C_S_3_tilde

# ╔═╡ 969e7fb0-99bd-4f6e-8be5-0dd04e354979
md"""
Since $S_2$ and $S_3$ are analytic in the unit disc they attain their maximum value on the boundary. It hence suffices to enclose the maximum for $|t| = 1$. Moreover, they are conjugate symmetric and it therefore suffices to bound them for $t = e^{i\theta}$ with $\theta \in [0, \pi]$.
"""

# ╔═╡ cb4549a2-a807-4b92-b2ff-6a6249bd2c77
@time S_2_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    abs(SRP.S(2, SRP.exppii(θ_div_π)))
end

# ╔═╡ ce427266-2654-41f0-9227-8f6fa54ea09f
@time S_3_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    rtol = 1e-3,
    maxevals = 5000,
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    abs(SRP.S(3, SRP.exppii(θ_div_π)))
end

# ╔═╡ 6f1af8dc-c7cf-45e1-976b-2b2106d4c5b8
md"""
We can now verify that the bound holds:
"""

# ╔═╡ 429393fe-79c2-41c2-9f03-5e1b60f7c194
S_2_bound <= Arb(C_S_2)

# ╔═╡ b845f335-a2ca-4d9a-b904-e054ef4b556e
S_3_bound <= Arb(C_S_3)

# ╔═╡ 3e288375-8972-45ad-b0fc-c76a841c4b1a
md"""
For the bounds after multiplication with $\frac{\rho^{1/2}}{\lambda^{1/2}}$ we directly enclose this factor as
"""

# ╔═╡ 66361fe4-88c4-4dcd-afac-fdaaeb947b80
sqrt_ρ_div_λ = let inv_N = Arb((0, 1 // N₀))
    SRP.c_N(inv_N) * sqrt(SRP.λ_app_div_λ(inv_N))
end

# ╔═╡ abdbe92c-7cad-4b8c-b50a-832fc616f865
md"""
Giving
"""

# ╔═╡ 7f3e1def-6986-43a5-b249-bfb866f2bf74
S_2_tilde_bound = sqrt_ρ_div_λ * S_2_bound

# ╔═╡ 07e837b0-b26c-44ac-b37c-051268853ef7
S_3_tilde_bound = sqrt_ρ_div_λ * S_3_bound

# ╔═╡ 1498b4a7-3a7f-446d-a4b4-45253dfd9ca9
md"""
Which we can verify satisfy the required bounds:
"""

# ╔═╡ 72ee9f48-7f79-4609-95c1-fff0f9f46aba
S_2_tilde_bound <= Arb(C_S_2_tilde)

# ╔═╡ 6a79fbd3-d7d7-4650-ad37-95b0343eeef7
S_3_tilde_bound <= Arb(C_S_3_tilde)

# ╔═╡ 624ba6af-d903-4722-8da4-46e17c7ebcdd
md"""
## Bounding $F_N$
We are in this case interested in showing that for $|t| \leq 1$ we have

$$\left|\frac{\rho^{1/2}}{\lambda^{1/2}}F_N(t)\right| \leq C_{F_N,0},$$

$$\left|\frac{\rho^{1/2}}{\lambda^{1/2}}F_N(t) - 1\right| \leq \frac{C_{F_N,2}}{N^2},$$

$$\left| \frac{\rho^{1/2}}{\lambda^{1/2}} \left(F_N(t) - \frac{S_2(t)}{N^2}\right) - 1 \right| \leq \frac{C_{F_N,3}}{N^3},$$

and

$$\left| \frac{\rho^{1/2}}{\lambda^{1/2}} \left(F_N(t) - \frac{S_2(t)}{N^2} - \frac{S_3(t)}{N^3} \right) - 1 \right| \leq \frac{C_{F_N,4}}{N^4}$$

With the constants given by
"""

# ╔═╡ 236a83c9-eba1-454d-ad8f-e044db582be6
C_F_N_0 = SRP.C_F_N_0

# ╔═╡ e74caacb-4be1-4ac1-a76b-7cc54ab7dc53
C_F_N_2 = SRP.C_F_N_2

# ╔═╡ e5801191-00ab-4325-9aa0-ce4953cdbf75
C_F_N_3 = SRP.C_F_N_3

# ╔═╡ 61843d44-005e-41f0-93e9-b8c07390e5f0
C_F_N_4 = SRP.C_F_N_4

# ╔═╡ 5f83199e-c80e-4a4c-a073-6acbbd72755f
md"""
We start by verifying the bound with $C_{F_N,4}$. Similar to before the function we are bounding is analytic on the unit disc and it therefore suffices to bound it on the boundary. Moreover, it is conjugate symmetric and it therefore suffices to bound it for $t = e^{i\theta}$ with $\theta \in [0, \pi]$.

To enclose

$$\frac{\rho^{1/2}}{\lambda^{1/2}} \left(F_N(t) - \frac{S_2(t)}{N^2} - \frac{S_3(t)}{N^3} \right) - 1$$

we compute a Taylor model of it of degree 3. By construction this Taylor model has the first four coefficients equal to zero, and the function is bounded by the remainder term divided by $N^4$.
"""

# ╔═╡ 1d8a5243-ea78-428d-9f89-0969a86d6c78
# This doesn't depend on t, so we can precompute it
sqrt_ρ_div_λ_model = let inv_N = Arb((0, 1 // N₀))
    ArbTaylorModel(inv_N, Arb(0), degree = 5) do inv_N
        SRP.c_N(inv_N) * sqrt(SRP.λ_app_div_λ(inv_N))
    end
end

# ╔═╡ c40e6da2-ac47-4a24-a420-ca8e673dfc25
function F_N_4(z)
    F_N_model = SRP.F_N_model(N₀, z)
    # Subtract N^-2 and N^-3 terms
    F_N_model.p[2] = 0
    F_N_model.p[3] = 0

    model = SRP.truncate(sqrt_ρ_div_λ_model * F_N_model - 1, degree = 3)

    # Verify that terms in the required expansion contain zero
    @assert all(Arblib.contains_zero.(model.p[0:3]))

    # Return the remainder term of the Taylor model
    model.p[end]
end

# ╔═╡ 8060aa4e-8689-4539-9009-4763902d11b7
@time F_N_4_bound = ArbExtras.maximum_enclosure(
    Arf(0),
    Arf(1),
    degree = -1,
    ubound_tol = Arb(C_F_N_4),
    abs_value = true,
    threaded = true,
    verbose = true,
) do θ_div_π
    abs(F_N_4(SRP.exppii(θ_div_π)))
end

# ╔═╡ 08097f8b-cce5-41e7-a311-008ff7aa4020
md"""
We can now verify the bound:
"""

# ╔═╡ 0839f91c-80a7-41f6-b239-5e79d26642d3
F_N_4_bound <= Arb(C_F_N_4)

# ╔═╡ bed41f2f-4e9c-4e37-8a56-34cfe6617c95
md"""
For the other cases we note that

$$\left|\frac{\rho^{1/2}}{\lambda^{1/2}}F_N(t)\right|
    \leq 1 + \left|\frac{\rho^{1/2}}{\lambda^{1/2}}\frac{S_2(t)}{N^3}\right| + \left|\frac{\rho^{1/2}}{\lambda^{1/2}}\frac{S_3(t)}{N^3}\right| + \frac{C_{F_N,4}}{N^4},
      \leq 1 + \frac{\tilde{C}_{S_{2}}}{N_{0}^{2}} + \frac{\tilde{C}_{S_{3}}}{N_{0}^{3}} + \frac{C_{F_N,4}}{N_{0}^{4}},$$

$$\left|\frac{\rho^{1/2}}{\lambda^{1/2}}F_N(t) - 1\right|
    \leq \left|\frac{\rho^{1/2}}{\lambda^{1/2}}\frac{S_2(t)}{N^3}\right| + \left|\frac{\rho^{1/2}}{\lambda^{1/2}}\frac{S_3(t)}{N^3}\right| + \frac{C_{F_N,4}}{N^4},
      \leq \left(\tilde{C}_{S_{2}} + \frac{\tilde{C}_{S_{3}}}{N_{0}} + \frac{C_{F_N,4}}{N_{0}^{2}}\right)\frac{1}{N^{2}}$$

and

$$\left| \frac{\rho^{1/2}}{\lambda^{1/2}} \left(F_N(t) - \frac{S_2(t)}{N^2} \right) - 1 \right|
    \leq \left|\frac{\rho^{1/2}}{\lambda^{1/2}}\frac{S_3(t)}{N^3}\right| + \frac{C_{F_N,4}}{N^4}
      \leq \left(\tilde{C}_{S_{3}} + \frac{C_{F_N,4}}{N_{0}}\right)\frac{1}{N^{3}}.$$

This gives us
"""

# ╔═╡ 78df064e-9b05-40ea-b56f-0653f08d1a02
F_N_0_bound = 1 + Arb(C_S_2_tilde) / N₀^2 + Arb(C_S_3_tilde) / N₀^3 + Arb(C_F_N_4) / N₀^4

# ╔═╡ 37ecc954-5f29-466f-832b-5018f55573fd
F_N_2_bound = Arb(C_S_2_tilde) + Arb(C_S_3_tilde) / N₀ + Arb(C_F_N_4) / N₀^2

# ╔═╡ a5633f97-395e-4d9d-98e9-9c1a5a5c19c8
F_N_3_bound = Arb(C_S_3_tilde) + Arb(C_F_N_4) / N₀

# ╔═╡ f68c3b8b-fada-445e-a9f1-85bc7ac659f4
md"""
We can now verify that these satisfy the required bounds:
"""

# ╔═╡ 4e37eb1b-0b60-4e38-a846-3604de33ec47
F_N_0_bound <= Arb(C_F_N_0)

# ╔═╡ 7e0502bc-cafd-41b3-a9ad-cb74c8b621d2
F_N_2_bound <= Arb(C_F_N_2)

# ╔═╡ bd593615-c270-4a5c-94af-920b4fc5494a
F_N_3_bound <= Arb(C_F_N_3)

# ╔═╡ 95276755-a29b-4365-8304-7e96b6f6c4b6
md"""
## Bounding $J_0^{(8)}$
TODO
"""

# ╔═╡ Cell order:
# ╟─0823e7e3-1b4d-409e-9668-150149d8aa1e
# ╠═9658d218-adae-11f0-1a53-675306a5a87f
# ╟─8c455525-a5ec-43d1-9afb-b31c8668a15e
# ╠═9362bb77-b3d5-463f-9962-229f832b3497
# ╟─f4ceb785-492f-4294-983e-c29360febd5f
# ╟─b037f618-1929-4c93-a452-5606cc59bb41
# ╠═abca81b1-bfdc-4b56-bae7-3d279ec14a29
# ╟─4711295e-0858-4833-b8e1-12384856ecad
# ╠═626e33ce-ece6-4066-af07-98bcb6c053c8
# ╟─a57c2ae1-85f3-4940-8d9b-2f4680fd79ef
# ╠═7c4055e6-60b9-4085-b9c8-7fe4da7c75d7
# ╟─a317e983-b641-4ca5-99a3-2b58fd0625c5
# ╠═c3c77e89-9059-45d6-a1e4-ab5f4a06d584
# ╟─c3264d17-3129-4cb7-a4d9-5e7acca2c2f8
# ╠═e6deace5-c025-482e-bc74-b8f66da238dd
# ╠═4b113a97-293f-4644-87fd-f536014caa71
# ╟─8a899965-2cce-4ce2-aca1-a9953a5ad9a8
# ╠═1f71ede2-30f8-4eb5-9faf-00137e2d64e8
# ╠═c97d4288-2112-48ed-889b-b00aaf4d925d
# ╟─969e7fb0-99bd-4f6e-8be5-0dd04e354979
# ╠═cb4549a2-a807-4b92-b2ff-6a6249bd2c77
# ╠═ce427266-2654-41f0-9227-8f6fa54ea09f
# ╟─6f1af8dc-c7cf-45e1-976b-2b2106d4c5b8
# ╠═429393fe-79c2-41c2-9f03-5e1b60f7c194
# ╠═b845f335-a2ca-4d9a-b904-e054ef4b556e
# ╟─3e288375-8972-45ad-b0fc-c76a841c4b1a
# ╠═66361fe4-88c4-4dcd-afac-fdaaeb947b80
# ╟─abdbe92c-7cad-4b8c-b50a-832fc616f865
# ╠═7f3e1def-6986-43a5-b249-bfb866f2bf74
# ╠═07e837b0-b26c-44ac-b37c-051268853ef7
# ╟─1498b4a7-3a7f-446d-a4b4-45253dfd9ca9
# ╠═72ee9f48-7f79-4609-95c1-fff0f9f46aba
# ╠═6a79fbd3-d7d7-4650-ad37-95b0343eeef7
# ╟─624ba6af-d903-4722-8da4-46e17c7ebcdd
# ╠═236a83c9-eba1-454d-ad8f-e044db582be6
# ╠═e74caacb-4be1-4ac1-a76b-7cc54ab7dc53
# ╠═e5801191-00ab-4325-9aa0-ce4953cdbf75
# ╠═61843d44-005e-41f0-93e9-b8c07390e5f0
# ╟─5f83199e-c80e-4a4c-a073-6acbbd72755f
# ╠═1d8a5243-ea78-428d-9f89-0969a86d6c78
# ╠═c40e6da2-ac47-4a24-a420-ca8e673dfc25
# ╠═8060aa4e-8689-4539-9009-4763902d11b7
# ╟─08097f8b-cce5-41e7-a311-008ff7aa4020
# ╠═0839f91c-80a7-41f6-b239-5e79d26642d3
# ╟─bed41f2f-4e9c-4e37-8a56-34cfe6617c95
# ╠═78df064e-9b05-40ea-b56f-0653f08d1a02
# ╠═37ecc954-5f29-466f-832b-5018f55573fd
# ╠═a5633f97-395e-4d9d-98e9-9c1a5a5c19c8
# ╟─f68c3b8b-fada-445e-a9f1-85bc7ac659f4
# ╠═4e37eb1b-0b60-4e38-a846-3604de33ec47
# ╠═7e0502bc-cafd-41b3-a9ad-cb74c8b621d2
# ╠═bd593615-c270-4a5c-94af-920b4fc5494a
# ╟─95276755-a29b-4365-8304-7e96b6f6c4b6
