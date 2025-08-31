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
    import ProgressLogging: @withprogress, @logprogress

    setprecision(Arb, 128)
end

# ╔═╡ c45a43db-337d-4ce7-afeb-61a9f50ba398
md"""
# Proof of Lemma 2.6
"""

# ╔═╡ b20ebf73-ffca-41f6-a2be-13d56881ad25
N₀ = 64

# ╔═╡ 4de76364-83fd-4889-91ff-cf0ae92ffdea
md"""
We want to prove that for $|z| = 1$ we have the following bounds:

$$|b_2(z)| \leq 3.5,\ |b_3(z)| \leq 2.5,\ |b_4(z)| \leq 10,\ |b_5(z)| \leq 12,$$

and for $N \geq$ $N₀ have

$$|T_2(z)| \leq 3.5,\ |T_4(z)| \leq 15,\ |T_6(z)| \leq 50.$$
"""

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
	SRP.b_2(SRP.exppii(θ_div_π))
end

# ╔═╡ e0ccdc24-1cd5-4967-86a4-226423d33fb5
b3s = tmap(θs_div_πs) do θ_div_π
	SRP.b_3(SRP.exppii(θ_div_π))
end

# ╔═╡ f8466db6-f0d4-44cc-b369-a587f8565bc7
b4s = tmap(θs_div_πs) do θ_div_π
	SRP.b_4(SRP.exppii(θ_div_π))
end

# ╔═╡ 95b752e2-5927-4235-babb-5a2198f172aa
b5s = tmap(θs_div_πs) do θ_div_π
	SRP.b_5(SRP.exppii(θ_div_π))
end

# ╔═╡ c25aea0b-f364-4c25-9fa8-8315fcc319cb
T6s = let T_6 = SRP.T_6_bound(N₀)
	tmap(θs_div_πs) do θ_div_π
		T_6(SRP.exppii(θ_div_π))
	end
end

# ╔═╡ 0209a6b1-f5ba-4ba6-9b5f-2f2cbdd2c643
let
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel = L"\theta",
		ylabel = L"b_2(e^{i\theta})",
	)
	scatterlines!(ax, θs, b2s)
	scatterlines!(ax, θs, abs.(b2s))
	fig
end

# ╔═╡ 1441945f-d537-4ac6-b21d-cfcc62c28cef
let
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel = L"\theta",
		ylabel = L"b_3(e^{i\theta})",
	)
	scatterlines!(ax, θs, b3s)
	scatterlines!(ax, θs, abs.(b3s))
	fig
end

# ╔═╡ 0c002f46-d981-4461-a2fe-0ceaf79ece35
let
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel = L"\theta",
		ylabel = L"b_4(e^{i\theta})",
	)
	scatterlines!(ax, θs, b4s)
	scatterlines!(ax, θs, abs.(b4s))
	fig
end

# ╔═╡ 0d186cc4-b8cd-4336-aa8d-e2ffdc37c80a
let
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel = L"\theta",
		ylabel = L"b_5(e^{i\theta})",
	)
	scatterlines!(ax, θs, b5s)
	scatterlines!(ax, θs, abs.(b5s))
	fig
end

# ╔═╡ 1f4ebea3-1b83-49d4-aabc-cc57b0bba1af
let
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel = L"\theta",
		ylabel = L"T_6(e^{i\theta})",
	)
	band!(ax, θs, lbound.(T6s), ubound.(T6s))
	#scatterlines!(ax, θs, lbound.(T6s))
	#scatterlines!(ax, θs, ubound.(T6s))
	fig
end

# ╔═╡ e6f37ebc-4884-4028-a937-93d71c27e3d0
md"""
## Proofs
"""

# ╔═╡ 0a8e3984-1047-4f89-8f16-d47b89fcd6cf
ArbExtras.maximum_enclosure(
	Arf(0),
	Arf(1),
	degree = -1,
	rtol = 1e-3,
	abs_value = true,
	threaded = true,
	verbose = true,
) do θ_div_π
	SRP.b_2(SRP.exppii(θ_div_π))
end

# ╔═╡ fd627d82-96e2-4f57-923b-dd879acf5bc6
ArbExtras.maximum_enclosure(
	Arf(0),
	Arf(1),
	degree = -1,
	rtol = 1e-3,
	abs_value = true,
	threaded = true,
	verbose = true,
) do θ_div_π
	SRP.b_3(SRP.exppii(θ_div_π))
end

# ╔═╡ b3f8fb15-2864-4479-8b2a-74586b2c8e8b
ArbExtras.maximum_enclosure(
	Arf(0),
	Arf(1),
	degree = -1,
	rtol = 1e-3,
	abs_value = true,
	threaded = true,
	verbose = true,
) do θ_div_π
	SRP.b_4(SRP.exppii(θ_div_π))
end

# ╔═╡ 00ac7282-9dac-4289-8d5a-5bf5dd539878
ArbExtras.maximum_enclosure(
	Arf(0),
	Arf(1),
	degree = -1,
	rtol = 1e-3,
	abs_value = true,
	threaded = true,
	verbose = true,
) do θ_div_π
	SRP.b_5(SRP.exppii(θ_div_π))
end

# ╔═╡ 763878e0-5a56-4c0e-985f-f0521cab30eb
let T_6 = SRP.T_6_bound(N₀)
	ArbExtras.maximum_enclosure(
		Arf(0),
		Arf(1),
		degree = -1,
		rtol = 1e-3,
		ubound_tol = 1.01ubound(abs(T_6(Acb(1)))),
		abs_value = true,
		threaded = true,
		verbose = true,
	) do θ_div_π
		T_6(SRP.exppii(θ_div_π))
	end
end

# ╔═╡ Cell order:
# ╟─c45a43db-337d-4ce7-afeb-61a9f50ba398
# ╠═cc7c7dea-83ab-11f0-04b8-a7eb5cc4403a
# ╠═b20ebf73-ffca-41f6-a2be-13d56881ad25
# ╠═4de76364-83fd-4889-91ff-cf0ae92ffdea
# ╟─e707c5eb-acd8-43fc-8fd1-1e0f9683d96b
# ╠═0aa4c803-768e-4b65-ba2c-c17f400d17d9
# ╠═28ef9722-2d9c-4ad3-ac78-baec417e9a45
# ╠═5d2dd24d-e9bd-434a-9a34-c4b9c130f1e6
# ╠═e0ccdc24-1cd5-4967-86a4-226423d33fb5
# ╠═f8466db6-f0d4-44cc-b369-a587f8565bc7
# ╠═95b752e2-5927-4235-babb-5a2198f172aa
# ╠═c25aea0b-f364-4c25-9fa8-8315fcc319cb
# ╠═0209a6b1-f5ba-4ba6-9b5f-2f2cbdd2c643
# ╟─1441945f-d537-4ac6-b21d-cfcc62c28cef
# ╟─0c002f46-d981-4461-a2fe-0ceaf79ece35
# ╠═0d186cc4-b8cd-4336-aa8d-e2ffdc37c80a
# ╠═1f4ebea3-1b83-49d4-aabc-cc57b0bba1af
# ╟─e6f37ebc-4884-4028-a937-93d71c27e3d0
# ╠═0a8e3984-1047-4f89-8f16-d47b89fcd6cf
# ╠═fd627d82-96e2-4f57-923b-dd879acf5bc6
# ╠═b3f8fb15-2864-4479-8b2a-74586b2c8e8b
# ╠═00ac7282-9dac-4289-8d5a-5bf5dd539878
# ╠═763878e0-5a56-4c0e-985f-f0521cab30eb
