### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 9384ba14-8538-11f0-0c3c-6f2e20c313cc
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

# ╔═╡ 93d81d7f-9ceb-4c5b-88e6-aa209faac94d
md"""
# Proof of Lemma 2.7
"""

# ╔═╡ 755359f0-227d-496f-b1bc-851e3579e149
md"""
## Goal
Let $g(w) = J_0(w\sqrt{\lambda})$, we want to prove that

$$|g'''(w)| \leq C_{g'''},$$

for $w \in [0, \infty)$, where $C_{g'''}$ is given by:
"""

# ╔═╡ 279cd351-92ec-40a7-889e-980444f43fed
C_gd3 = SRP.C_gd3

# ╔═╡ 72009879-ad52-46db-931c-97af0db81870
md"""
## Implementation
Note that

$$g'''(w) = -\frac{1}{4}\lambda^{3/2}\left(3J_1(w\sqrt{\lambda}) - J_3(w\sqrt{\lambda})\right).$$

We can implement this as:
"""

# ╔═╡ f756cd4a-a496-4366-b47b-71f14532431b
λ = SRP.λ_disc()

# ╔═╡ 26942198-2078-4d35-8444-f63d6e9e94aa
g_d3(w) = -λ^(3 // 2) * (3besselj1(w * sqrt(λ) - besselj(Arb(3), w * sqrt(λ)))) / 4

# ╔═╡ d0364f27-acbf-41d7-a744-2e653db06914
md"""
## Proof
We split the interval $[0, \infty)$ into two parts, $[0, 5]$ and $(5, \infty)$. On the interval $[0, 5]$ we directly enclose the maximum:
"""

# ╔═╡ dd258019-77b9-4da2-874e-896162b18c78
g_d3_bound_0_5 = ArbExtras.maximum_enclosure(
    g_d3,
    Arf(0),
    Arf(5),
    rtol = 1e-5,
    abs_value = true,
    verbose = true,
)

# ╔═╡ 4b00cb82-1e2c-448c-b832-b94ed22d92c6
md"""
On the interval $(5, \infty)$ we use the bound

$$|J_\nu(x)| \leq 0.7858x^{-1/3},$$

which holds for $x > 0$. From this we immediately get

$$|g'''(w)| \leq 0.7858\lambda^{4/3}w^{-1/3} \leq 0.7858\lambda^{4/3}5^{-1/3}$$

Where this last value is given by:
"""

# ╔═╡ 8e11a8f8-3b11-4b6d-a8ad-c6a368e20abb
g_d3_bound_5_inf = Arb("0.7858") * λ^(4 // 3) * Arb(5)^(-1 // 3)

# ╔═╡ 9001f1d9-8977-4fbf-8f7d-1091ab04ffe9
md"""
For inclusing in the paper we print it using fewer digits.
"""

# ╔═╡ f50b7428-db8a-4511-a8c1-b8bdd8eb5dec
string(g_d3_bound_5_inf, digits = 5)

# ╔═╡ e62a9ae7-d492-4a03-9c9d-9736ce9de35a
md"""
It is then straight forward to verify that the bounds on $[0, 5]$ and $(5, \infty)$ satisfy the required inequality:
"""

# ╔═╡ bc230116-438e-4b7f-874e-1b18a83bab92
@assert_proof g_d3_bound_0_5 <= Arb(C_gd3)

# ╔═╡ da0960f3-4793-47a8-a1cf-a813b55e2a9e
@assert_proof g_d3_bound_5_inf <= Arb(C_gd3)

# ╔═╡ 9dabbf48-5575-40f2-93cb-35be265f423a
md"""
## Plots
We can visualize the above bounds in the plot below. Note that this plot is not part of the proof, it is only intended to give a visualization of the results.
"""

# ╔═╡ c1cf75a7-2e6d-43a4-8715-009c872be9ef
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"w")
    colormap = :tab10
    colorrange = (1, 10)
    lines!(
        ax,
        range(0, 20, 1000),
        g_d3 ∘ Arb,
        label = L"g'''(w)",
        color = 1;
        colormap,
        colorrange,
    )
    lines!(
        ax,
        range(0, 20, 1000),
        w -> 0.7858 * λ^(4 // 3) * w^(-1 // 3),
        label = L"\pm 0.7858\lambda^{4/3}w^{-1/3}",
        color = 2;
        colormap,
        colorrange,
    )
    lines!(ax, range(0, 20, 1000), w -> -0.7858 * λ^(4 // 3) * w^(-1 // 3), color = :orange)
    hlines!(
        ax,
        [-Arb(C_gd3), Arb(C_gd3)],
        label = L"\pm C_g",
        color = 3;
        colormap,
        colorrange,
    )
    axislegend(ax, position = :rb)
    fig
end

# ╔═╡ Cell order:
# ╟─93d81d7f-9ceb-4c5b-88e6-aa209faac94d
# ╠═9384ba14-8538-11f0-0c3c-6f2e20c313cc
# ╟─755359f0-227d-496f-b1bc-851e3579e149
# ╠═279cd351-92ec-40a7-889e-980444f43fed
# ╟─72009879-ad52-46db-931c-97af0db81870
# ╠═f756cd4a-a496-4366-b47b-71f14532431b
# ╠═26942198-2078-4d35-8444-f63d6e9e94aa
# ╟─d0364f27-acbf-41d7-a744-2e653db06914
# ╠═dd258019-77b9-4da2-874e-896162b18c78
# ╟─4b00cb82-1e2c-448c-b832-b94ed22d92c6
# ╠═8e11a8f8-3b11-4b6d-a8ad-c6a368e20abb
# ╟─9001f1d9-8977-4fbf-8f7d-1091ab04ffe9
# ╠═f50b7428-db8a-4511-a8c1-b8bdd8eb5dec
# ╟─e62a9ae7-d492-4a03-9c9d-9736ce9de35a
# ╠═bc230116-438e-4b7f-874e-1b18a83bab92
# ╠═da0960f3-4793-47a8-a1cf-a813b55e2a9e
# ╟─9dabbf48-5575-40f2-93cb-35be265f423a
# ╟─c1cf75a7-2e6d-43a4-8715-009c872be9ef
