### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 4c62a864-0a2d-4968-afad-9a3d3b7a9113
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

# ╔═╡ 755affb2-aeb8-11f1-2002-c15dda69fa11
md"""
# Proof of Lemma 3.2
"""

# ╔═╡ 57230f56-0f3a-423c-8546-eada78657387
md"""
## Goal
We want to prove that for $N \geq N_0 + 1$ we have

$$\lambda_N < \lambda_2(\mathcal{P}_N),$$

where $N_0$ is given by:
"""

# ╔═╡ f6dc8c30-2269-455f-ac1b-58ff5b0ff401
N₀ = SRP.N₀

# ╔═╡ 3847253f-f089-4afa-9c8f-d3e4a802d37e
md"""
## Proof

The proof is separated into two parts. The first is to check $N=3$ and $N=4$ with the known explicit expressions for their eigenvalues. For the remaining cases $5 \leq N \leq N_0+1$, we want to prove that

$$\lambda_N < \lambda_2(\mathbb{D}_{C_5}),$$

where $\lambda_2(\mathbb{D}_{C_5})$ denotes the second eigenvalue of $\mathbb{D}_{C_5}$ and $C_5$ is the circumradius of $\mathcal{P}_5$.

### Step 1- Check condition for $N=3,4$

For the equilateral triangle of area $\pi$,

$$\lambda_1(\mathcal P_3)=\frac{4\pi}{\sqrt3},
\qquad
\lambda_2(\mathcal P_3)=\frac{28\pi}{3\sqrt3}.$$

For the square of area $\pi$,

$$\lambda_1(\mathcal P_4)=2\pi,
\qquad
\lambda_2(\mathcal P_4)=5\pi.$$
"""

# ╔═╡ 289f041b-727c-4f79-af4e-26825fafbfae
λ₁_triangle = 4Arb(π) / sqrt(Arb(3))

# ╔═╡ b75e5ed8-a7b6-40f1-a4a8-be2fe789d919
λ₂_triangle = 28Arb(π) / (3 * sqrt(Arb(3)))

# ╔═╡ f68ee0ea-2b57-4942-8a0e-e8661e76c105
λ₁_square = 2Arb(π)

# ╔═╡ 33144ac2-4543-4650-ba06-c56d046f222e
λ₂_square = 5Arb(π)

# ╔═╡ e31dd340-228b-41c8-a00a-aceacbd98cfe
@assert_proof λ₁_triangle < λ₂_triangle

# ╔═╡ 45f71b49-f759-4c6e-b07c-64abee5e114b
@assert_proof λ₁_square < λ₂_square

# ╔═╡ a1a0b048-7386-4f59-b2a5-4930032604ee
md"""
### Step 2 - Compute approximations for $5\leq N \leq N_0+1$

We compute approximations of $\lambda_{1}(\mathcal{P}_N)$ as well as the associated eigenfunction for $5 \leq N \leq N_0 + 1$.

For $N \geq 12$ we use precomputed approximations from [the
  repository](https://github.com/David-Berghaus/master-thesis-data) associated with[Computation of Laplacian eigenvalues of two-dimensional shapes with dihedral symmetry](http://dx.doi.org/10.1007/s10444-024-10138-3), for $5 \leq N \leq 11$ we compute an approximation on the fly.
"""

# ╔═╡ e4f8270f-9ac9-4a4b-9e68-794f4c4418f9
Ns = 5:1:(N₀+1)

# ╔═╡ 19d41475-6127-4c5c-8058-de4ee9c87d4d
domains = RegularPolygon{Arb}.(Ns)

# ╔═╡ 2513f9b0-8863-4d6d-9e62-665264cd5b47
us, λs_approx = let
    us = Eigenfunction.(domains)

    M = 4 # Number of terms to use in the approximations

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
        if u.domain.N <= 11
            # For these values of N using Float64 instead of BigFloat for the
            # computations seems to be more stable
            SRP.sigma!(u, λ_approx, M, qr_eltype = Float64)
        else
            SRP.sigma!(u, λ_approx, M)
        end
    end

    us, λs_approx
end

# ╔═╡ d1d86bca-7f87-410f-83d3-0483e940d8b8
md"""
Next we compute rigorous enclosures of an eigenvalue for $5 \leq N \leq N_0 + 1$.
"""

# ╔═╡ e396fad2-76c2-4d3c-ba1a-0427c75a48af
@time λs = let
    progress = Threads.Atomic{Int}(0)
    @withprogress tmap(Arb, us, λs_approx, scheduler = :greedy) do u, λ_approx
        λ = SRP.eigenvalue_enclosure(u, λ_approx, threaded = false)
        Threads.atomic_add!(progress, 1)
        @logprogress progress[] / length(Ns)
        λ
    end
end

# ╔═╡ ed1a52cd-440e-4eb5-89f0-e8fc7bb1b260
md"""
We check that all enclosures were successfully computed.
"""

# ╔═╡ ce461b32-9c52-4cb4-a6ee-529ac876eb0c
@assert_proof all(isfinite, λs)

# ╔═╡ b3990cc6-618e-4457-9710-140338e418df
md"""
### Step 3 - Compute $\lambda_2(\mathbb{D}_{C_5})$
We now compute an enclosure of $\lambda_2(\mathbb{D}_{C_5})$. By scaling we get that

$$\lambda_2(\mathbb{D}_{C_5}) = \frac{\lambda_2(\mathbb{D})}{C_5^2}.$$

To compute $\lambda_2(\mathbb{D})$ we use that it is given by

$$\lambda_2(\mathbb{D}) = j_{1,1}^2,$$

where $j_{1,1}$ denotes the first positive zero of $J_1$.
"""

# ╔═╡ 92cccb5d-f2e9-4021-888a-961e3cf18100
md"""
First we compute $C_5$ from the circumradius formula of a regular $N$-polygon of area $\pi$:

$$C_N = \sqrt{\frac{\pi}{\frac{N}{2} \sin\left(\frac{2\pi}{N}\right)}}.$$
"""

# ╔═╡ 13326883-4f99-4ead-a176-51aeaf2a2806
function circumradius_area_pi(N::Integer)
    θ = 2Arb(π) / Arb(N)
    sqrt(2Arb(π) / (Arb(N) * sin(θ)))
end

# ╔═╡ ba2b2c52-e0fc-4513-852d-e9568b4d3dba
C_5 = circumradius_area_pi(5)

# ╔═╡ 8dec50c4-e906-4419-aa3e-62a306341a11
md"""
To compute $j_{1,1}$ we first verify that $J_1$ is strictly increasing on the interval $[0, 1]$, this ensures that the only zero on that interval is the one at zero.
"""

# ╔═╡ a4f11ff2-3c96-47a5-8681-969b8bedf523
Arblib.ispositive(ArbExtras.derivative_function(besselj1)(Arb((0, 1))))

# ╔═╡ 75aa04ae-8662-47a0-a46b-1081ed1ca04f
md"""
This proves that $j_{1,1} > 1$. Next we isolate all the roots on the interval $[1, 5]$:
"""

# ╔═╡ fe8662ec-d7a8-49b5-bf7f-716e45ddecf0
roots, flags = ArbExtras.isolate_roots(besselj1, Arf(1), Arf(5), depth = 20)

# ╔═╡ 992c68c7-95d1-479c-8d94-e08911eba7e4
md"""
We verify that it only found one root and that it was proved to be unique:
"""

# ╔═╡ c22247bd-775c-42c9-8401-8c8c5aa7102b
length(flags) == 1 && flags[1]

# ╔═╡ 09fd387f-b019-4535-9866-c5aa9658980e
md"""
We refine the enclosure of the root, giving an enclosure for $j_{1,1}$:
"""

# ╔═╡ 95e368db-53bc-4117-9a75-2f13f0b9c4f0
j_1_1 = ArbExtras.refine_root(besselj1, Arb(roots[1]))

# ╔═╡ defe50e3-1f63-4e15-849f-9435a21d71a4
md"""
We square it, to get an enclosure of the eigenvalue for $\mathbb{D}$:
"""

# ╔═╡ 2a64d7d8-d861-47eb-98fe-af28554714e6
λ₂ = j_1_1^2

# ╔═╡ 7d91f783-5977-49de-bb4f-56cc5031af57
md"""
Finally we scale the result to get an enclosure of $\lambda_2(\mathbb{D}_{C_5})$:
"""

# ╔═╡ 071a86a0-e5f0-4c44-88db-a2a802bffb34
λ₂_C_5 = λ₂ / Arb(C_5)^2

# ╔═╡ 3f857190-e15d-49e7-8659-7c48822b842b
md"""
We print a version with fewer digits for inclusion in the paper.
"""

# ╔═╡ ae014df2-f806-444d-81ad-61f5449fbcf5
string(λ₂_C_5, digits = 5)

# ╔═╡ 7e1f9d00-995f-4417-b131-54e9b6980b9a
md"""
### Step 4 - Verify inequality
Finally, we just verify the inequality:
"""

# ╔═╡ a3955c31-f647-441a-bc36-b4a07a124eaa
@assert_proof all(eachindex(λs)) do i
    λs[i] < λ₂_C_5
end

# ╔═╡ Cell order:
# ╟─755affb2-aeb8-11f1-2002-c15dda69fa11
# ╠═4c62a864-0a2d-4968-afad-9a3d3b7a9113
# ╟─57230f56-0f3a-423c-8546-eada78657387
# ╠═f6dc8c30-2269-455f-ac1b-58ff5b0ff401
# ╟─3847253f-f089-4afa-9c8f-d3e4a802d37e
# ╠═289f041b-727c-4f79-af4e-26825fafbfae
# ╠═b75e5ed8-a7b6-40f1-a4a8-be2fe789d919
# ╠═f68ee0ea-2b57-4942-8a0e-e8661e76c105
# ╠═33144ac2-4543-4650-ba06-c56d046f222e
# ╠═e31dd340-228b-41c8-a00a-aceacbd98cfe
# ╠═45f71b49-f759-4c6e-b07c-64abee5e114b
# ╟─a1a0b048-7386-4f59-b2a5-4930032604ee
# ╠═e4f8270f-9ac9-4a4b-9e68-794f4c4418f9
# ╠═19d41475-6127-4c5c-8058-de4ee9c87d4d
# ╠═2513f9b0-8863-4d6d-9e62-665264cd5b47
# ╟─d1d86bca-7f87-410f-83d3-0483e940d8b8
# ╠═e396fad2-76c2-4d3c-ba1a-0427c75a48af
# ╟─ed1a52cd-440e-4eb5-89f0-e8fc7bb1b260
# ╠═ce461b32-9c52-4cb4-a6ee-529ac876eb0c
# ╟─b3990cc6-618e-4457-9710-140338e418df
# ╟─92cccb5d-f2e9-4021-888a-961e3cf18100
# ╠═13326883-4f99-4ead-a176-51aeaf2a2806
# ╠═ba2b2c52-e0fc-4513-852d-e9568b4d3dba
# ╟─8dec50c4-e906-4419-aa3e-62a306341a11
# ╠═a4f11ff2-3c96-47a5-8681-969b8bedf523
# ╟─75aa04ae-8662-47a0-a46b-1081ed1ca04f
# ╠═fe8662ec-d7a8-49b5-bf7f-716e45ddecf0
# ╟─992c68c7-95d1-479c-8d94-e08911eba7e4
# ╠═c22247bd-775c-42c9-8401-8c8c5aa7102b
# ╟─09fd387f-b019-4535-9866-c5aa9658980e
# ╠═95e368db-53bc-4117-9a75-2f13f0b9c4f0
# ╟─defe50e3-1f63-4e15-849f-9435a21d71a4
# ╠═2a64d7d8-d861-47eb-98fe-af28554714e6
# ╟─7d91f783-5977-49de-bb4f-56cc5031af57
# ╠═071a86a0-e5f0-4c44-88db-a2a802bffb34
# ╟─3f857190-e15d-49e7-8659-7c48822b842b
# ╠═ae014df2-f806-444d-81ad-61f5449fbcf5
# ╟─7e1f9d00-995f-4417-b131-54e9b6980b9a
# ╠═a3955c31-f647-441a-bc36-b4a07a124eaa
