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

# ╔═╡ 7c2b770e-a59b-4f0a-af8a-fcf6c4c55324
N₀ = 64

# ╔═╡ bdc6dcff-50a0-49d7-98a2-00b24f5ecfe6
md"""
We want to prove that for $|z| = 1$ we have the following bounds:

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} d_1(z, t)V_4(t)\ dt\right| \leq 15,$$

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} d_2(z, t)V_3(t)\ dt\right| \leq 10,$$

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} d_3(z, t)V_2(t)\ dt\right| \leq 35$$

and for $N \geq$ $N₀ have

$$\left|\operatorname{Re}\int_0^z \frac{1}{t} K_4(z, t)V_1(t)\ dt\right| \leq 40.$$
"""

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

# ╔═╡ 7bc1266a-4f01-439b-81fb-5b1e86612ef2
let
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel = L"\theta",
		#ylabel = L"T_6(e^{i\theta})",
	)
	band!(ax, θs, lbound.(real(I_d_1_V_4)), ubound.(real(I_d_1_V_4)))
	hlines!(ax, [-15, 15])
	fig
end

# ╔═╡ 6187163b-832d-4c9f-862f-a3483505fda6
let
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel = L"\theta",
		#ylabel = L"T_6(e^{i\theta})",
	)
	band!(ax, θs, lbound.(real(I_d_2_V_3)), ubound.(real(I_d_2_V_3)))
	hlines!(ax, [-10, 10])
	fig
end

# ╔═╡ a7509228-24ea-4395-8bdc-57457c63e28b
let
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel = L"\theta",
		#ylabel = L"T_6(e^{i\theta})",
	)
	band!(ax, θs, lbound.(real(I_d_3_V_2)), ubound.(real(I_d_3_V_2)))
	hlines!(ax, [-35, 35])
	fig
end

# ╔═╡ 93cf97f7-54c1-4101-93b8-09890db8de86


# ╔═╡ Cell order:
# ╟─18d792d6-8c0c-4214-bd98-670611da3e1e
# ╠═fca89648-853b-11f0-3116-5b610b900eb4
# ╠═7c2b770e-a59b-4f0a-af8a-fcf6c4c55324
# ╠═bdc6dcff-50a0-49d7-98a2-00b24f5ecfe6
# ╟─b6096ef5-2ce3-4519-97db-6ccb566037f0
# ╠═ef2ecb67-904b-4488-8e47-b5e2fc46eaab
# ╠═d24e7344-a2cd-4ba3-a5cf-1290fe285cfc
# ╠═cf4e04c7-26f2-400b-a885-c3e495e42f78
# ╠═67e1ad91-7da5-4985-b8c3-ca1cc84a47fb
# ╠═a820d24b-9387-4ef5-a959-a4baaf93fdb2
# ╠═7bc1266a-4f01-439b-81fb-5b1e86612ef2
# ╠═6187163b-832d-4c9f-862f-a3483505fda6
# ╠═a7509228-24ea-4395-8bdc-57457c63e28b
# ╠═93cf97f7-54c1-4101-93b8-09890db8de86
