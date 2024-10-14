### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(
                Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                "AbstractPlutoDingetjes",
            )].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6fae0fd4-3232-11ef-1c63-21b4630e06dc
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using SpectralRegularPolygons
    using GLMakie
    using Arblib
    using ArbExtras
    using OhMyThreads
    using PlutoUI
    using SpecialFunctions

    import SpectralRegularPolygons as SRP
    import ProgressLogging: @withprogress, @logprogress

    setprecision(BigFloat, 128)
    setprecision(Arb, 128)
end

# ╔═╡ a3f86580-1615-4109-98ac-a0eb403c9a78
md"""
# Monotonicity of eigenvalue
This notebook contains the computer assisted part of Proposition 3.1 in the accompanying paper. More precisely it proves that for $N_{0} = 26$ we have

$$\lambda_{1}(\mathcal{P}_3) > \lambda_{1}(\mathcal{P}_4) >
\cdots > \lambda_{1}(\mathcal{P}_{N_{0} - 1}) > \lambda_{1}(\mathcal{P}_{N_{0}}),$$

where $\lambda_{1}(\mathcal{P}_N)$ denotes the first eigenvalue of the regular $N$-gon with area $\pi$. Furthermore, if 

$$q_N := \frac{\lambda_1(\mathcal{P}_N)}{\lambda_1(\mathcal{P}_{N+1})}$$

it asserts that

$$q_{3} > q_{4} > \dots > q_N > q_{N+1} > \dots >  q_{N_{0} - 1} > q_{N_{0}}.$$

Note that $\lambda_{1}(\mathcal{P}_3) ) = \frac{4\pi}{\sqrt{3}}$ and $\lambda_{1}(\mathcal{P}_4) ) = 2\pi$ are both known exactly. There is therefore no need to compute enclosures of these eigenvalues using the Method of Particular Solutions, we do however still have to verify that they satisfy the required inequalities.
"""

# ╔═╡ 04daab37-a7ea-472a-950a-da2c8a45ce71
N₀ = 26

# ╔═╡ 4377e7e0-d63c-4be2-9dbe-70f86515509d
md"""
## Construct approximations
"""

# ╔═╡ 674a8185-fe50-4f70-a0a5-e4cbe95460c4
Ns = 5:1:N₀+1

# ╔═╡ 5a88affb-12d4-46ed-a3dc-c02e00d1ecdf
domains = RegularPolygon{Arb}.(Ns)

# ╔═╡ 2a08cf06-8258-4efb-9aae-3bc66c617c37
us, λs_approx = let
    us = Eigenfunction.(domains)

    M = 3

    λs_approx = tmap(us, SRP.get_eigenvalue_approximation.(Arb, Ns)) do u, λ
        if !isfinite(λ)
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

# ╔═╡ b35d0bfa-2f2b-43e9-b59a-e94e0e23e2cb
@time λs = let
    progress = Threads.Atomic{Int}(0)
    @withprogress tmap(Arb, us, λs_approx, scheduler = :greedy) do u, λ_approx
        if u.domain.N == 3
            λ = 4Arb(π) / sqrt(Arb(3))
        elseif u.domain.N == 4
            λ = 2Arb(π)
        else
            λ = SRP.eigenvalue_enclosure(u, λ_approx, threaded = false)
        end
        @logprogress Threads.atomic_add!(progress, 1)[] / length(Ns)
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
The code above computes the eigenvalues for $N \geq 5$. We also want to check the conditions for $N = 3$ and $N = 4$, so we create a vector with those two added.
"""

# ╔═╡ 468dccc1-36f7-4441-a410-7d5c99341d63
λs_full = [4Arb(π) / sqrt(Arb(3)); 2Arb(π); λs]

# ╔═╡ e44b149a-6c13-4464-a39c-17317e9a4ff5
md"""
Next we verify that the inequality $\lambda_{1}(\mathcal{P}_{N}) > \lambda_{1}(\mathcal{P}_{N + 1})$ holds.
"""

# ╔═╡ 827a5af8-ac00-414c-b79b-e8f98a907969
all(eachindex(λs_full)[1:end-1]) do i
    λs_full[i] > λs_full[i+1]
end

# ╔═╡ 99bd9807-92c7-4723-be12-a66700e04587
md"""
Finally we compute the values for $q_N$, and verify that $q_{N} > q_{N + 1}$ holds.
"""

# ╔═╡ a65680cb-6cc5-445d-9658-4cf94c79306a
qs_full = λs_full[1:end-1] ./ λs_full[2:end]

# ╔═╡ bcbd6354-5212-4efc-a42c-0f70a1e49bc2
all(eachindex(qs_full)[1:end-1]) do i
    qs_full[i] > qs_full[i+1]
end

# ╔═╡ d4143768-2eb4-4d52-9321-7daf30823c90
md"""
## Produce plots
"""

# ╔═╡ 8b6b7d10-467d-4dcd-a800-2bbe0e9573a0
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save_figures CheckBox(default = false))
"""

# ╔═╡ e771c0f7-9e80-4c92-ab37-fc8f3e9866df
let
    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1, 1], xlabel = L"N", ylabel = L"\lambda")

    scatterlines!(ax, Ns, λs)

    save_figures && save("figures/eigenvalues.png", fig, px_per_unit = 4)

    fig
end

# ╔═╡ 1190007b-9b1a-4063-ba96-f98b93137b10
let
    a1s = map(u -> u.vertex_expansion.coefficients[1], us)
    a2s = map(u -> u.vertex_expansion.coefficients[2], us)

    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1, 1], xlabel = L"N")

    scatterlines!(ax, Ns, a1s, label = L"a_1")
    scatterlines!(ax, Ns, a2s, label = L"a_2")
    axislegend(ax)

    save_figures && save("figures/coefficients.png", fig, px_per_unit = 4)

    fig
end

# ╔═╡ bda98e53-cb8d-4a8d-a3c9-f521dae9a187
let
    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1, 1], xlabel = L"N", yscale = log10)

    λs_diff = λs[1:end-1] - λs[2:end]
    radius_sum = Arblib.radius.(λs[1:end-1]) + Arblib.radius.(λs[2:end])

    scatterlines!(
        ax,
        Ns[1:end-1],
        λs_diff,
        label = L"\lambda_{1}(\mathcal{P}_N) - \lambda_{1}(\mathcal{P}_{N + 1})",
    )
    scatterlines!(ax, Ns[1:end-1], radius_sum, label = "Error bound")
    axislegend(ax)

    save_figures && save("figures/difference.png", fig, px_per_unit = 4)

    fig
end

# ╔═╡ ea97c15b-b03f-4b6b-8e80-73ff12eaa215
md"""
### Figures not appearing the in the paper
"""

# ╔═╡ 8e96f457-c878-441e-b70b-c47980439963
let
    λ_circle = ArbExtras.refine_root(besselj0, Arb((sqrt(Arf(5.7)), sqrt(Arf(5.8)))))^2

    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1, 1], xlabel = L"N", yscale = log10)

    scatterlines!(
        ax,
        Ns,
        λs_approx .- λ_circle,
        label = L"\lambda_{1}(\mathcal{P}_N) - \lambda_{1}(\mathbb{D})",
    )
    scatterlines!(ax, Ns, Arblib.radius.(λs), label = "Error bound")
    axislegend(ax)

    fig
end

# ╔═╡ Cell order:
# ╟─a3f86580-1615-4109-98ac-a0eb403c9a78
# ╠═6fae0fd4-3232-11ef-1c63-21b4630e06dc
# ╠═04daab37-a7ea-472a-950a-da2c8a45ce71
# ╟─4377e7e0-d63c-4be2-9dbe-70f86515509d
# ╠═674a8185-fe50-4f70-a0a5-e4cbe95460c4
# ╠═5a88affb-12d4-46ed-a3dc-c02e00d1ecdf
# ╠═2a08cf06-8258-4efb-9aae-3bc66c617c37
# ╟─c2150b3f-00d7-4400-9cd5-a12f21c9e469
# ╠═b35d0bfa-2f2b-43e9-b59a-e94e0e23e2cb
# ╟─540eabb1-4dd7-4e80-890a-5577bb64cf36
# ╠═d5b9ea36-fee2-4b8b-a090-16c2ecf308f3
# ╟─fac49c09-f365-4fe7-a05a-7b6bfa306e50
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
# ╟─bda98e53-cb8d-4a8d-a3c9-f521dae9a187
# ╟─ea97c15b-b03f-4b6b-8e80-73ff12eaa215
# ╟─8e96f457-c878-441e-b70b-c47980439963
