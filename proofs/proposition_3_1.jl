### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 6fae0fd4-3232-11ef-1c63-21b4630e06dc
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
    import SpectralRegularPolygons: RegularPolygon, Eigenfunction
    import ProgressLogging: @withprogress, @logprogress

    setprecision(BigFloat, 128)
    setprecision(Arb, 128)
end

# ╔═╡ 0a9e2d34-1a0f-4ce9-9ed8-1097f1c2700c
md"""
# Proof of Proposition 3.1
"""

# ╔═╡ a3f86580-1615-4109-98ac-a0eb403c9a78
md"""
This notebook contains the computer assisted part of Proposition 3.1. More precisely it proves that for $N_{0}$ given by
"""

# ╔═╡ 2ee8254b-2a81-4047-84d5-bea626ad6b2c
N₀ = SRP.N₀

# ╔═╡ cd3086d7-ac55-4ed2-834d-856314972790
md"""
we have

$$\lambda_{1}(\mathcal{P}_3) > \lambda_{1}(\mathcal{P}_4) >
\cdots > \lambda_{1}(\mathcal{P}_{N_{0} - 1}) > \lambda_{1}(\mathcal{P}_{N_{0}}),$$

where $\lambda_{1}(\mathcal{P}_N)$ denotes the first eigenvalue of the regular $N$-gon with area $\pi$. Furthermore, if

$$q_N = \frac{\lambda_1(\mathcal{P}_N)}{\lambda_1(\mathcal{P}_{N+1})}$$

it asserts that

$$q_{3} > q_{4} > \dots > q_N > q_{N+1} > \dots >  q_{N_{0} - 1} > q_{N_{0}}.$$

Note that $\lambda_{1}(\mathcal{P}_3) ) = \frac{4\pi}{\sqrt{3}}$ and $\lambda_{1}(\mathcal{P}_4) ) = 2\pi$ are both known exactly. There is therefore no need to compute enclosures of these eigenvalues using the Method of Particular Solutions, we do however still have to verify that they satisfy the required inequalities.
"""

# ╔═╡ 4377e7e0-d63c-4be2-9dbe-70f86515509d
md"""
## Construct approximations
We start by computing approximations of $\lambda_{1}(\mathcal{P}_N)$ as well as the associated eigenfunction for $5 \leq N \leq N_0 + 1$.

For $N \geq 12$ we use precomputed approximations, for $5 \leq N \leq 11$ we compute an approximation on the fly.
"""

# ╔═╡ 674a8185-fe50-4f70-a0a5-e4cbe95460c4
Ns = 5:1:(N₀+1)

# ╔═╡ 5a88affb-12d4-46ed-a3dc-c02e00d1ecdf
domains = RegularPolygon{Arb}.(Ns)

# ╔═╡ 2a08cf06-8258-4efb-9aae-3bc66c617c37
us, λs_approx = let
    us = Eigenfunction.(domains)

    M = 4

    λs_approx = tmap(us, SRP.get_eigenvalue_approximation.(Arb, Ns)) do u, λ
        if !isfinite(λ)
            # Comptue approximation on the fly
            @assert 5 <= u.domain.N <= 11
            λ_upper = SRP.get_eigenvalue_approximation(Arb, 4)
            λ_lower = SRP.get_eigenvalue_approximation(Arb, 12)
            SRP.mps(u, λ_lower, λ_upper, M, qr_eltype = Float64)
        else
            λ
        end
    end

    tforeach(us, λs_approx, scheduler = :greedy) do u, λ_approx
        if u.domain.N <= 12
            # For these values of N using Float64 instead of BigFloat for the 
            # computations seems to be more stable
            SRP.sigma!(u, λ_approx, M, qr_eltype = Float64)
        else
            SRP.sigma!(u, λ_approx, M)
        end
    end

    us, λs_approx
end

# ╔═╡ c2150b3f-00d7-4400-9cd5-a12f21c9e469
md"""
## Compute enclosures
"""

# ╔═╡ b8b858a8-b6f0-4780-8d38-6fce97c28d43
md"""
Next we compute rigorous enclosures of $\lambda_{1}(\mathcal{P}_N)$ for $5 \leq N \leq N_0 + 1$. 
"""

# ╔═╡ b35d0bfa-2f2b-43e9-b59a-e94e0e23e2cb
@time λs = let
    progress = Threads.Atomic{Int}(0)
    @withprogress tmap(Arb, us, λs_approx, scheduler = :greedy) do u, λ_approx
        λ = SRP.eigenvalue_enclosure(u, λ_approx, threaded = false)
        Threads.atomic_add!(progress, 1)
        @logprogress progress[] / length(Ns)
        λ
    end
end

# ╔═╡ 540eabb1-4dd7-4e80-890a-5577bb64cf36
md"""
## Verify proposition
With the enclosures computed, the next step is to check that they satisfy the required conditions. First we check that all enclosures were succesfully computed.
"""

# ╔═╡ d5b9ea36-fee2-4b8b-a090-16c2ecf308f3
all(isfinite, λs)

# ╔═╡ fac49c09-f365-4fe7-a05a-7b6bfa306e50
md"""
The code above computeted the eigenvalues for $5 \leq N \leq N_0 + 1$. We also want to check the conditions for $N = 3$ and $N = 4$, so we create a vector with those two added.
"""

# ╔═╡ c115ab1e-cd48-4c95-aebf-3c04715c28fa
Ns_full = 3:Ns[end]

# ╔═╡ 468dccc1-36f7-4441-a410-7d5c99341d63
λs_full = [4Arb(π) / sqrt(Arb(3)); 2Arb(π); λs]

# ╔═╡ e44b149a-6c13-4464-a39c-17317e9a4ff5
md"""
Next we verify that the inequality $\lambda_{1}(\mathcal{P}_{N}) > \lambda_{1}(\mathcal{P}_{N + 1})$ holds.
"""

# ╔═╡ 827a5af8-ac00-414c-b79b-e8f98a907969
all(eachindex(λs_full)[1:(end-1)]) do i
    λs_full[i] > λs_full[i+1]
end

# ╔═╡ 99bd9807-92c7-4723-be12-a66700e04587
md"""
Finally we compute the values for $q_N$ for $3 \leq N \leq N_0$, and verify that $q_{N} > q_{N + 1}$ holds.
"""

# ╔═╡ a65680cb-6cc5-445d-9658-4cf94c79306a
qs_full = λs_full[1:(end-1)] ./ λs_full[2:end]

# ╔═╡ bcbd6354-5212-4efc-a42c-0f70a1e49bc2
all(eachindex(qs_full)[1:(end-1)]) do i
    qs_full[i] > qs_full[i+1]
end

# ╔═╡ d4143768-2eb4-4d52-9321-7daf30823c90
md"""
## Produce plots
Finally we produce some plots for inclusion in the paper. See the paper for more details about what exactly they show.
"""

# ╔═╡ 8b6b7d10-467d-4dcd-a800-2bbe0e9573a0
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save_figures CheckBox(default = false))
"""

# ╔═╡ e771c0f7-9e80-4c92-ab37-fc8f3e9866df
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N", ylabel = L"\lambda")

    scatterlines!(ax, Ns, λs)

    save_figures && save("figures/eigenvalues-small-N.pdf", fig)

    fig
end

# ╔═╡ 1190007b-9b1a-4063-ba96-f98b93137b10
let
    a2s = map(u -> u.interior_expansion.coefficients[2], us)
    a2s_scaling = besselj.(Arb.(Ns), Arb(1))

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N", yscale = log10)

    scatterlines!(ax, Ns, abs.(a2s) .* a2s_scaling, label = L"|a_2| J_N(1)")
    axislegend(ax)

    save_figures && save("figures/coefficient-a2-small-N.pdf", fig)

    fig
end

# ╔═╡ ea84defb-8199-41c9-a427-486eb67c881a
let
    b1s = map(u -> u.vertex_expansion.coefficients[1], us)
    b2s = map(u -> u.vertex_expansion.coefficients[2], us)

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N")

    scatterlines!(ax, Ns, b1s, label = L"b_1", marker = :circle)
    scatterlines!(ax, Ns, b2s, label = L"b_2", marker = :cross)
    axislegend(ax)

    save_figures && save("figures/coefficients-bs-small-N.pdf", fig)

    fig
end

# ╔═╡ bda98e53-cb8d-4a8d-a3c9-f521dae9a187
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N", yscale = log10)

    λs_diff = λs_full[1:(end-1)] - λs_full[2:end]
    radius_sum = radius.(λs_full[1:(end-1)]) + radius.(λs_full[2:end])

    scatterlines!(
        ax,
        Ns_full[1:(end-2)],
        λs_diff[1:(end-1)],
        label = L"\lambda_{1}(\mathcal{P}_N) - \lambda_{1}(\mathcal{P}_{N + 1})",
    )
    scatterlines!(
        ax,
        Ns_full[2:(end-2)],
        radius_sum[2:(end-1)],
        label = "Error bound",
        marker = :cross,
    )
    axislegend(ax)

    save_figures && save("figures/difference-eigenvalues-small-N.pdf", fig)

    fig
end

# ╔═╡ f594e504-9527-4c43-8e15-96713c15f9de
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N", yscale = log10)

    qs_diff = qs_full[1:(end-1)] - qs_full[2:end]
    radius_sum = radius.(qs_full[1:(end-1)]) + radius.(qs_full[2:end])

    scatterlines!(ax, Ns_full[1:(end-2)], qs_diff, label = L"q_{N} - q_{N + 1}")
    scatterlines!(
        ax,
        Ns_full[2:(end-2)],
        radius_sum[2:end],
        label = "Error bound",
        marker = :cross,
    )
    axislegend(ax)

    save_figures && save("figures/difference-q-N-small-N.pdf", fig)

    fig
end

# ╔═╡ ea97c15b-b03f-4b6b-8e80-73ff12eaa215
md"""
### Figures not appearing the in the paper
"""

# ╔═╡ 01fce4da-8495-4de0-809a-272171c66357
md"""
This figure shows the distance between the eigenvalue of the polygon and the disc, as well as the radius of the enclosure for the eigenvalue.
"""

# ╔═╡ 8e96f457-c878-441e-b70b-c47980439963
let
    λ_circle = ArbExtras.refine_root(besselj0, Arb((sqrt(Arf(5.7)), sqrt(Arf(5.8)))))^2

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"N", yscale = log10)

    scatterlines!(
        ax,
        Ns_full,
        λs_full .- λ_circle,
        label = L"\lambda_{1}(\mathcal{P}_N) - \lambda_{1}(\mathbb{D})",
    )
    scatterlines!(ax, Ns, Arblib.radius.(λs), label = "Error bound")
    axislegend(ax)

    fig
end

# ╔═╡ 511e9aa5-2322-4f6c-8278-f8db03a213ec
md"""
This figure shows the behavior of the approximate eigenfunction on the boundary of the polygon. Note that we only plot it along half of one boundary segment due to symmetry.
"""

# ╔═╡ 5045160f-3ac3-4774-98db-3cb93e4812de
let
    N = 5
    i = N - 4
    domain = domains[i]
    u = us[i]
    λ_approx = λs_approx[i]

    fig = Figure()
    ax = Axis(fig[1, 1], title = L"u_1(\mathcal{P}_%$N)")

    ts = range(Arb(0), 1, 100)[2:(end-1)]
    points = SpectralRegularPolygons.boundary_parameterized_symmetry.(Ref(domain), ts)
    values = u.(points, λ_approx)

    scatterlines!(ax, ts, values)

    fig
end

# ╔═╡ Cell order:
# ╟─0a9e2d34-1a0f-4ce9-9ed8-1097f1c2700c
# ╠═6fae0fd4-3232-11ef-1c63-21b4630e06dc
# ╟─a3f86580-1615-4109-98ac-a0eb403c9a78
# ╠═2ee8254b-2a81-4047-84d5-bea626ad6b2c
# ╟─cd3086d7-ac55-4ed2-834d-856314972790
# ╟─4377e7e0-d63c-4be2-9dbe-70f86515509d
# ╠═674a8185-fe50-4f70-a0a5-e4cbe95460c4
# ╠═5a88affb-12d4-46ed-a3dc-c02e00d1ecdf
# ╠═2a08cf06-8258-4efb-9aae-3bc66c617c37
# ╟─c2150b3f-00d7-4400-9cd5-a12f21c9e469
# ╟─b8b858a8-b6f0-4780-8d38-6fce97c28d43
# ╠═b35d0bfa-2f2b-43e9-b59a-e94e0e23e2cb
# ╟─540eabb1-4dd7-4e80-890a-5577bb64cf36
# ╠═d5b9ea36-fee2-4b8b-a090-16c2ecf308f3
# ╟─fac49c09-f365-4fe7-a05a-7b6bfa306e50
# ╠═c115ab1e-cd48-4c95-aebf-3c04715c28fa
# ╠═468dccc1-36f7-4441-a410-7d5c99341d63
# ╟─e44b149a-6c13-4464-a39c-17317e9a4ff5
# ╠═827a5af8-ac00-414c-b79b-e8f98a907969
# ╟─99bd9807-92c7-4723-be12-a66700e04587
# ╠═a65680cb-6cc5-445d-9658-4cf94c79306a
# ╠═bcbd6354-5212-4efc-a42c-0f70a1e49bc2
# ╟─d4143768-2eb4-4d52-9321-7daf30823c90
# ╟─8b6b7d10-467d-4dcd-a800-2bbe0e9573a0
# ╟─e771c0f7-9e80-4c92-ab37-fc8f3e9866df
# ╟─1190007b-9b1a-4063-ba96-f98b93137b10
# ╟─ea84defb-8199-41c9-a427-486eb67c881a
# ╟─bda98e53-cb8d-4a8d-a3c9-f521dae9a187
# ╟─f594e504-9527-4c43-8e15-96713c15f9de
# ╟─ea97c15b-b03f-4b6b-8e80-73ff12eaa215
# ╟─01fce4da-8495-4de0-809a-272171c66357
# ╟─8e96f457-c878-441e-b70b-c47980439963
# ╟─511e9aa5-2322-4f6c-8278-f8db03a213ec
# ╟─5045160f-3ac3-4774-98db-3cb93e4812de
