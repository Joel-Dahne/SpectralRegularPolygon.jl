### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ a8ef131c-85bb-11f0-28cf-7b4c670ed13b
begin
    using Pkg
    Pkg.activate("..", io = devnull)
    using SpectralRegularPolygons
    using Arblib
    using ArbExtras

    import SpectralRegularPolygons as SRP

    setprecision(Arb, 128)
end

# ╔═╡ bfb90022-f09a-46a0-b0aa-5d1ffdd7a908
md"""
# Proof of Lemma 2.22
"""

# ╔═╡ d2d34a83-67a8-4d71-8b06-171fefb55487
md"""
## Goal
We want to prove that for $N \geq N_0$ we have

$$\frac{\lambda_{\text{app}}(N)}{1 + \hat{\varepsilon}(N)} > \frac{\lambda_{\text{app}}(N+1)}{1 - \hat{\varepsilon}(N+1)}$$

and

$$\frac{\lambda_{\text{app}}(N)(1-\hat{\varepsilon}(N+1))}{(1+\hat{\varepsilon}(N))\lambda_{\text{app}}(N+1)} > \frac{\lambda_{\text{app}}(N+1)(1+\hat{\varepsilon}(N+2))}{(1-\hat{\varepsilon}(N+1))\lambda_{\text{app}}(N+2)},$$

where $N_0$ is given by:
"""

# ╔═╡ 1e814a20-d22d-4c0d-9003-6712a0a5d205
N₀ = SRP.N₀

# ╔═╡ 992bfff9-2fc8-4bf7-a3d7-3e5de4861881
md"""
## Proof
Let $P_1(N)$ and $P_2(N)$ denote the difference between the left and right hand side of the two inequalities. Our goal is then to show that $P_1(N) > 0$ and $P_2(N) > 0$ for $N \geq N_0$.

For the computations it is easier to work with $\nu = 1 / N$ instead of $N$. For that reason we introduce the functions $p_1$ and $p_2$, which we define by $p_1(\nu) = P_1(1 / \nu)$ and $p_2(\nu) = P_2(1 / \nu)$.

**Note:** In the package we mostly use the variable name `inv_N` to refer to `1 / N`, in this notebook we however use `ν` to more closely follow the notation of the proof in the paper.
"""

# ╔═╡ 1259b964-c46b-49ed-93d4-828df69a6677
md"""
With the above notation we have

$$p_1(\nu) = \frac{\lambda_{\text{app}}(1/\nu)}{1 + \hat{\varepsilon}(1/\nu)} - \frac{\lambda_{\text{app}}(1/\nu+1)}{1 - \hat{\varepsilon}(1/\nu+1)}$$

from which we factor out $\lambda$ and use `SRP.λ_app_div_λ` to enclose $\lambda_\text{app}/ \lambda$. Note that $(1 / \nu + 1)^{-1} = \nu / (1 + \nu)$
"""

# ╔═╡ 5dc8f540-22ee-4040-b91b-4810ca27627f
function p_1(ν)
    inv_inv_ν_p1 = ν / (1 + ν) # Enclosure of inv(inv(ν) + 1)

    return SRP.λ_disc() * (
        SRP.λ_app_div_λ(ν) / (1 + SRP.epsilon_hat(ν)) -
        SRP.λ_app_div_λ(inv_inv_ν_p1) / (1 - SRP.epsilon_hat(inv_inv_ν_p1))
    )
end

# ╔═╡ a38c6989-a3d7-4b72-9e97-aa4d53040dab
md"""
Furthermore we have

$$p_2(\nu) =\frac{\lambda_{\text{app}}(1/\nu)(1-\hat{\varepsilon}(1/\nu+1))}{(1+\hat{\varepsilon}(1/\nu))\lambda_{\text{app}}(1/\nu+1)} > \frac{\lambda_{\text{app}}(1/\nu+1)(1+\hat{\varepsilon}(1\nu+2))}{(1-\hat{\varepsilon}(1/\nu+1))\lambda_{\text{app}}(1/\nu+2)}$$

from which explicitly cancel a factor $\lambda$ between the numerator and denominator and use `SRP.λ_app_div_λ` to enclose $\lambda_\text{app}/ \lambda$. Note that $(1 / \nu + 1)^{-1} = \nu / (1 + \nu)$ and $(1 / \nu + 2)^{-1} = \nu / (1 + 2\nu)$
"""

# ╔═╡ 24b69dd3-b353-4ed1-9ef6-ed79f8a9ad3d
function p_2(ν)
    inv_inv_ν_p1 = ν / (1 + ν) # Enclosure of inv(inv(ν) + 1)
    inv_inv_ν_p2 = ν / (1 + 2ν) # Enclosure of inv(inv(ν) + 2)

    # Note that we use λ_app_div_λ since the λs cancel
    return SRP.λ_app_div_λ(ν) * (1 - SRP.epsilon_hat(inv_inv_ν_p1)) /
           ((1 + SRP.epsilon_hat(ν)) * SRP.λ_app_div_λ(inv_inv_ν_p1)) -
           SRP.λ_app_div_λ(inv_inv_ν_p1) * (1 + SRP.epsilon_hat(inv_inv_ν_p2)) /
           ((1 - SRP.epsilon_hat(inv_inv_ν_p1)) * SRP.λ_app_div_λ(inv_inv_ν_p2))
end

# ╔═╡ ecdfb51c-bd88-4d2b-86fd-d6c46340bf28
md"""
### Part 1
Let us start by proving that $p_1(\nu) > 0$.

From the proof in the paper we have that the first four terms in the expansion at $\nu = 0$ vanish. It therefore suffices to show that the fourth derivative of $p_1$ is positive on $[0, 1 / N_0]$ to ensure that $p_1$ is also positive.
"""

# ╔═╡ 8a2a6a79-ef3c-4efc-a75b-e3871c88a875
md"""
As a double check for that the first four terms in the expansion vanish we can verify that the expansion computed by the computer agrees with this. Computing the expansion we get:
"""

# ╔═╡ 76af3933-de21-45f0-9310-028eaf31465f
p_1_expansion = p_1(ArbSeries((0, 1), degree = 4))

# ╔═╡ bb882554-bd78-4ab2-acfc-e7fb702e0024
md"""
The first three terms in the expansion are computed to be exactly zero:
"""

# ╔═╡ b20f3e75-b475-4cac-b69d-b6611ac82363
iszero(p_1_expansion[0]) && iszero(p_1_expansion[1]) && iszero(p_1_expansion[2])

# ╔═╡ bfcf8109-a214-4460-b41a-35a714096f05
md"""
For the fourth term the computations cannot directly prove that it cancels, but we can verify that the enclosure contains zero:
"""

# ╔═╡ 0d216fba-d53a-47ba-886f-dc8902b0b0c2
Arblib.contains_zero(p_1_expansion[3])

# ╔═╡ e03b99aa-2802-470e-adb3-dbf6a74da074
md"""
What remains is to prove that the fourth derivative is positive on the entire interval $[0, 1 / N_0]$. We therefore enclose the minimum value of $p^{(4)}$ on this interval:
"""

# ╔═╡ f7c11696-e359-48c6-bf29-b799b3a26957
p_1_d4_minimum = ArbExtras.minimum_enclosure(
    ArbExtras.derivative_function(p_1, 4),
    Arf(0),
    ubound(Arb(1 // N₀)),
    verbose = true,
)

# ╔═╡ 232ee031-d05d-44f7-8ce9-fca770b67428
string(p_1_d4_minimum, digits = 5) # For inclusion in the paper

# ╔═╡ b2b04733-b22f-4283-b657-013bd02864ac
md"""
Finally we verify that the minimum is positive:
"""

# ╔═╡ 1b5615d5-ff94-41c6-8123-d11bdae8af08
@assert_proof Arblib.ispositive(p_1_d4_minimum)

# ╔═╡ 8184b4d2-7235-4cfe-8a34-a267ae170be7
md"""
### Part 2
"""

# ╔═╡ 4371e67c-6c2e-461c-bbb5-afcb8ad7f142
md"""
Next we prove that $p_2(\nu) > 0$

In this case we have from the proof in the paper that the first five terms in the expansion at $\nu = 0$ vanish. However, contrary to the situation for $p_1$, $p_2^{(5)}$ is not positive on the entire interval $[0, 1 / N_0]$. In particular it is negative at $1 / N_0$:
"""

# ╔═╡ 75b5d8f6-3f16-4ce0-965f-a9b8cee80196
string(ArbExtras.derivative_function(p_2, 5)(Arb(1 // N₀)), digits = 10)

# ╔═╡ a4811176-133b-4559-84b7-022e473df33f
md"""
We therefore split $[0, 1 / N_0]$ into $[0, a]$ and $[a, 1 / N_0]$, with:
"""

# ╔═╡ dbf0f411-344f-4d6c-b031-b5c61022cbd4
a = Arf(1 // 1024)

# ╔═╡ 0b7faaac-134f-4aa0-b880-b16de83083cf
md"""
On $[0, a]$ we want to prove that the fifth derivative is positive. We therefore enclose the minimum value of $p_2^{(5)}$ on this interval:
"""

# ╔═╡ a14c2a9a-4561-4420-b7aa-fe386bc84ea2
p_2_d5_minimum_0_a = ArbExtras.minimum_enclosure(
    ArbExtras.derivative_function(p_2, 5),
    Arf(0),
    a,
    verbose = true,
)

# ╔═╡ 47efff6f-f2cf-4ba7-8892-36948ec2f645
string(p_2_d5_minimum_0_a, digits = 5) # For inclusion in the paper

# ╔═╡ b8208fd1-664a-49fd-9f62-6f132941e010
md"""
We can verify that the minimum is positive:
"""

# ╔═╡ 8b5219a8-f395-49f4-a279-84fd0f744c92
@assert_proof Arblib.ispositive(p_2_d5_minimum_0_a)

# ╔═╡ ccaf3416-57dc-4c31-825f-104c408faa79
md"""
On the interval $[a, 1 / N_0]$ we directly enclose the minimum value of $p_2$ and verify that it is positive.
"""

# ╔═╡ 94b731fb-019c-43e1-b0b3-519d0f801529
p_2_minimum_a_inv_N₀ =
    ArbExtras.minimum_enclosure(p_2, a, ubound(Arb(1 // N₀)), verbose = true)

# ╔═╡ b90dddf4-2498-4015-9579-866d079a831b
string(p_2_minimum_a_inv_N₀, digits = 5) # For inclusion in the paper

# ╔═╡ cdeb73d6-600e-49b4-8dd9-2329d4b9dc9d
@assert_proof Arblib.ispositive(p_2_minimum_a_inv_N₀)

# ╔═╡ 2f58afa7-1597-4486-a14b-09128da921ff
md"""
This concludes the proof!

Similar to for $p_1$ we can also double check that the computations agree with the first five terms in the expansion at $\nu = 0$ vanishing.
"""

# ╔═╡ c427a061-0bf9-4e40-947e-2bd093ad686a
p_2_expansion = p_2(ArbSeries((0, 1), degree = 5))

# ╔═╡ dced2ca8-b61c-4d2f-a3cb-f0302376fe09
md"""
The first three terms in the expansion are computed to be exactly zero:
"""

# ╔═╡ 8374aeda-b25a-4a9a-b392-47b423892b86
iszero(p_2_expansion[0]) && iszero(p_2_expansion[1]) && iszero(p_2_expansion[2])

# ╔═╡ dd0ceaef-7f67-4d0e-ac35-b88ea25d430c
md"""
For the fourth and fifth term the computations cannot directly prove that they cancel, but we can verify that the enclosures contain zero:
"""

# ╔═╡ 218d06e7-c755-41b3-bfaa-1dfd6a3ccded
Arblib.contains_zero(p_2_expansion[3]) && Arblib.contains_zero(p_2_expansion[4])

# ╔═╡ Cell order:
# ╟─bfb90022-f09a-46a0-b0aa-5d1ffdd7a908
# ╠═a8ef131c-85bb-11f0-28cf-7b4c670ed13b
# ╟─d2d34a83-67a8-4d71-8b06-171fefb55487
# ╠═1e814a20-d22d-4c0d-9003-6712a0a5d205
# ╟─992bfff9-2fc8-4bf7-a3d7-3e5de4861881
# ╟─1259b964-c46b-49ed-93d4-828df69a6677
# ╠═5dc8f540-22ee-4040-b91b-4810ca27627f
# ╟─a38c6989-a3d7-4b72-9e97-aa4d53040dab
# ╠═24b69dd3-b353-4ed1-9ef6-ed79f8a9ad3d
# ╟─ecdfb51c-bd88-4d2b-86fd-d6c46340bf28
# ╟─8a2a6a79-ef3c-4efc-a75b-e3871c88a875
# ╠═76af3933-de21-45f0-9310-028eaf31465f
# ╟─bb882554-bd78-4ab2-acfc-e7fb702e0024
# ╠═b20f3e75-b475-4cac-b69d-b6611ac82363
# ╟─bfcf8109-a214-4460-b41a-35a714096f05
# ╠═0d216fba-d53a-47ba-886f-dc8902b0b0c2
# ╟─e03b99aa-2802-470e-adb3-dbf6a74da074
# ╠═f7c11696-e359-48c6-bf29-b799b3a26957
# ╠═232ee031-d05d-44f7-8ce9-fca770b67428
# ╟─b2b04733-b22f-4283-b657-013bd02864ac
# ╠═1b5615d5-ff94-41c6-8123-d11bdae8af08
# ╟─8184b4d2-7235-4cfe-8a34-a267ae170be7
# ╟─4371e67c-6c2e-461c-bbb5-afcb8ad7f142
# ╠═75b5d8f6-3f16-4ce0-965f-a9b8cee80196
# ╟─a4811176-133b-4559-84b7-022e473df33f
# ╠═dbf0f411-344f-4d6c-b031-b5c61022cbd4
# ╟─0b7faaac-134f-4aa0-b880-b16de83083cf
# ╠═a14c2a9a-4561-4420-b7aa-fe386bc84ea2
# ╠═47efff6f-f2cf-4ba7-8892-36948ec2f645
# ╟─b8208fd1-664a-49fd-9f62-6f132941e010
# ╠═8b5219a8-f395-49f4-a279-84fd0f744c92
# ╟─ccaf3416-57dc-4c31-825f-104c408faa79
# ╠═94b731fb-019c-43e1-b0b3-519d0f801529
# ╠═b90dddf4-2498-4015-9579-866d079a831b
# ╠═cdeb73d6-600e-49b4-8dd9-2329d4b9dc9d
# ╟─2f58afa7-1597-4486-a14b-09128da921ff
# ╠═c427a061-0bf9-4e40-947e-2bd093ad686a
# ╟─dced2ca8-b61c-4d2f-a3cb-f0302376fe09
# ╠═8374aeda-b25a-4a9a-b392-47b423892b86
# ╟─dd0ceaef-7f67-4d0e-ac35-b88ea25d430c
# ╠═218d06e7-c755-41b3-bfaa-1dfd6a3ccded
