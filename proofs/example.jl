### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 94a83dc6-320f-11ef-0fc1-47f8a7c333f2
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using SpectralRegularPolygons
    using GLMakie
    using Arblib
    using GeometryBasics
    using OhMyThreads

    import SpectralRegularPolygons as SRP

    setprecision(BigFloat, 128)
    setprecision(Arb, 128)
end

# ╔═╡ 8a477d5e-15ec-45af-a7bb-f9a9d2edd319
md"""
## Set parameters and construct domain
"""

# ╔═╡ 8e404533-aaea-48f3-87ac-074c2c2af210
N = 12

# ╔═╡ 9177dcec-02e6-4087-81ee-a91d2a091aa0
T = Arb

# ╔═╡ faba5d43-9ebc-47a2-9a5d-8da107cc4b47
domain = RegularPolygon{T}(N)

# ╔═╡ b840780a-f0cc-465b-a8ce-c4458bf7f4eb
md"""
## Get approximate eigenvalue and eigenfunction
"""

# ╔═╡ 98f95f01-3fe1-4d41-8602-57d0adc114a5
_λ = SRP.get_eigenvalue_approximation(T, N)

# ╔═╡ 7d66ba8b-f78b-42b3-a31a-fb69b0c81f8b
u, λ = let
    u = Eigenfunction(domain)

    M = 8

    if !isfinite(_λ)
        @assert 5 <= N <= 11
        λ_upper = SRP.get_eigenvalue_approximation(T, 4)
        λ_lower = SRP.get_eigenvalue_approximation(T, 12)
        λ = SRP.mps(u, λ_lower, λ_upper, M, qr_eltype = Float64)
    else
        λ = _λ
    end

    if N <= 12
        SRP.sigma!(u, λ, M, qr_eltype = Float64)
    else
        SRP.sigma!(u, λ, M)
    end
    u, λ
end

# ╔═╡ 6744e813-f571-4ed1-856e-943b8e87bbd2
md"""
## Enclose eigenvalue
"""

# ╔═╡ 0294e116-8d75-4d52-b6c2-8cd7a4282c7b
λ_enclosure = SRP.eigenvalue_enclosure(u, λ, verbose = true)

# ╔═╡ 77f1398c-523c-4b18-9c1a-6f357600c860
md"""
## Plotting
"""

# ╔═╡ bf9ca063-0f37-4b27-9cf2-714614e59010
SRP.plot_domain(domain)

# ╔═╡ 9db897cc-f282-4054-8f93-07891db9f470
SRP.plot_eigenfunction(u, λ, num_grid_points = 50)

# ╔═╡ 8c2dee50-13e4-4ebf-ab02-c4d66c27ddf3
SRP.plot_boundary(u, λ)

# ╔═╡ Cell order:
# ╠═94a83dc6-320f-11ef-0fc1-47f8a7c333f2
# ╟─8a477d5e-15ec-45af-a7bb-f9a9d2edd319
# ╠═8e404533-aaea-48f3-87ac-074c2c2af210
# ╠═9177dcec-02e6-4087-81ee-a91d2a091aa0
# ╠═faba5d43-9ebc-47a2-9a5d-8da107cc4b47
# ╟─b840780a-f0cc-465b-a8ce-c4458bf7f4eb
# ╠═98f95f01-3fe1-4d41-8602-57d0adc114a5
# ╠═7d66ba8b-f78b-42b3-a31a-fb69b0c81f8b
# ╟─6744e813-f571-4ed1-856e-943b8e87bbd2
# ╠═0294e116-8d75-4d52-b6c2-8cd7a4282c7b
# ╟─77f1398c-523c-4b18-9c1a-6f357600c860
# ╠═bf9ca063-0f37-4b27-9cf2-714614e59010
# ╠═9db897cc-f282-4054-8f93-07891db9f470
# ╠═8c2dee50-13e4-4ebf-ab02-c4d66c27ddf3
