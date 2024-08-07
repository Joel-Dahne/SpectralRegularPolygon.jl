### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 6fae0fd4-3232-11ef-1c63-21b4630e06dc
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using SpectralRegularPolygons
    using GLMakie
    using Arblib
    using ArbExtras
    using OhMyThreads
    using SpecialFunctions

    import SpectralRegularPolygons as SRP
    import ProgressLogging: @withprogress, @logprogress

    setprecision(BigFloat, 128)
    setprecision(Arb, 128)
end

# ╔═╡ 8f82b94e-ff18-430e-a669-ed9de7a5d97c
T = Arb

# ╔═╡ 674a8185-fe50-4f70-a0a5-e4cbe95460c4
Ns = 5:1:24

# ╔═╡ 5a88affb-12d4-46ed-a3dc-c02e00d1ecdf
domains = RegularPolygon{T}.(Ns)

# ╔═╡ be9a629d-2fb0-4e93-8b6d-ccd0c2e23a3f
_λs = SRP.get_eigenvalue_approximation.(T, Ns)

# ╔═╡ 2a08cf06-8258-4efb-9aae-3bc66c617c37
us, λs = let
    us = Eigenfunction.(domains)

    M = 3

    λs = tmap(us, _λs) do u, _λ
        if !isfinite(_λ)
            @assert 5 <= u.domain.N <= 11
            λ_upper = SRP.get_eigenvalue_approximation(T, 4)
            λ_lower = SRP.get_eigenvalue_approximation(T, 12)
            SRP.mps(u, λ_lower, λ_upper, M, qr_eltype = Float64)
        else
            _λ
        end
    end

    tforeach(us, λs, scheduler = :greedy) do u, λ
        if u.domain.N <= 12
            SRP.sigma!(u, λ, M, qr_eltype = Float64)
        else
            SRP.sigma!(u, λ, M)
        end
    end

    us, λs
end

# ╔═╡ b35d0bfa-2f2b-43e9-b59a-e94e0e23e2cb
@time λs_enclosures = let
    progress = Threads.Atomic{Int}(0)
    @withprogress tmap(Arb, us, λs, scheduler = :greedy) do u, λ
        λ_enclosure = SRP.eigenvalue_enclosure(u, λ, threaded = false)
        Threads.atomic_add!(progress, 1)
        @logprogress progress[] / length(Ns)
        λ_enclosure
    end
end

# ╔═╡ d5b9ea36-fee2-4b8b-a090-16c2ecf308f3
all(isfinite, λs_enclosures)

# ╔═╡ ea943df1-b84f-4e61-99a8-96a26df1bca2
all(λs_enclosures[1:end-1] .> λs_enclosures[2:end])

# ╔═╡ e771c0f7-9e80-4c92-ab37-fc8f3e9866df
plot(Ns, λs_enclosures)

# ╔═╡ 8e96f457-c878-441e-b70b-c47980439963
let
    λ_circle = ArbExtras.refine_root(besselj0, Arb((sqrt(Arf(5.7)), sqrt(Arf(5.8)))))^2
    fig = GLMakie.Figure()
    axis = GLMakie.Axis(fig[1, 1], yscale = log10)
    scatterlines!(axis, Ns, λs .- λ_circle)
    scatterlines!(axis, Ns, Arblib.radius.(λs_enclosures))
    fig
end

# ╔═╡ 194bc594-4cf3-4a87-acba-5b0b59f820aa
let
    fig = GLMakie.Figure()
    axis = GLMakie.Axis(fig[1, 1], yscale = log10)

    scatter!(axis, Ns, Arblib.radius.(λs_enclosures))

    fig
end

# ╔═╡ bda98e53-cb8d-4a8d-a3c9-f521dae9a187
let
    fig = GLMakie.Figure()
    axis = GLMakie.Axis(fig[1, 1], yscale = log10)

    λs_diff = λs[1:end-1] - λs[2:end]
    radius_sum =
        Arblib.radius.(λs_enclosures[1:end-1]) + Arblib.radius.(λs_enclosures[2:end])

    scatter!(axis, Ns[1:end-1], λs_diff)
    scatter!(axis, Ns[1:end-1], radius_sum)

    fig
end

# ╔═╡ Cell order:
# ╠═6fae0fd4-3232-11ef-1c63-21b4630e06dc
# ╠═8f82b94e-ff18-430e-a669-ed9de7a5d97c
# ╠═674a8185-fe50-4f70-a0a5-e4cbe95460c4
# ╠═5a88affb-12d4-46ed-a3dc-c02e00d1ecdf
# ╠═be9a629d-2fb0-4e93-8b6d-ccd0c2e23a3f
# ╠═2a08cf06-8258-4efb-9aae-3bc66c617c37
# ╠═b35d0bfa-2f2b-43e9-b59a-e94e0e23e2cb
# ╠═d5b9ea36-fee2-4b8b-a090-16c2ecf308f3
# ╠═ea943df1-b84f-4e61-99a8-96a26df1bca2
# ╠═e771c0f7-9e80-4c92-ab37-fc8f3e9866df
# ╠═8e96f457-c878-441e-b70b-c47980439963
# ╠═194bc594-4cf3-4a87-acba-5b0b59f820aa
# ╠═bda98e53-cb8d-4a8d-a3c9-f521dae9a187
