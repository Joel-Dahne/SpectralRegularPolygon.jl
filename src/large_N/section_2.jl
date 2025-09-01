"""
    λ_disc()

Compute an enclosure of the first eigenvalue of the disc.

**PROVE:** We need to verify that this indeed is the first zero. We
could do this with ArbExtras.isolate_roots, but it requires more
computations if we do this every time.
"""
λ_disc() = ArbExtras.refine_root(besselj0, Arb((sqrt(Arf(5.7)), sqrt(Arf(5.8)))))^2

"""
    λ_approx(inv_N::Union{Arb,ArbSeries})

Compute an enclosure of
```
λ * (1 + 4zeta(3) / N^3 + (12 - 2λ) * zeta(5) / N^5)
```
Note that this takes as input `inv(N)` and not `N`.
"""
function λ_approx(inv_N::Union{Arb,ArbSeries})
    λ = λ_disc()
    return λ * (1 + 4zeta(Arb(3)) * inv_N^3 + (12 - 2λ) * zeta(Arb(5)) * inv_N^5)
end

"""
    λ_approx_div_λ(inv_N::Union{Arb,ArbSeries})

Compute an enclosure of [`λ_approx`](@ref) divided by `λ`, i.e.
```
1 + 4zeta(3) / N^3 + (12 - 2λ) * zeta(5) / N^5
```
Note that this takes as input `inv(N)` and not `N`.
"""
λ_approx_div_λ(inv_N::Union{Arb,ArbSeries}) =
    1 + 4zeta(Arb(3)) * inv_N^3 + (12 - 2λ_disc()) * zeta(Arb(5)) * inv_N^5

"""
    sqrt_λ_approx_div_λ_remainder_N3(N₀::Int)

Compute an enclosure of
```
(sqrt(λ_approx / λ) - 1) * N^3
```
that is valid for all `N >= N₀`.
"""
function sqrt_λ_approx_div_λ_remainder_N3(N₀::Int)
    N_max = 100

    # Compute enclosure for N₀ <= N <= N_max - 1
    λ = λ_disc()
    values = map(N₀:(N_max-1)) do N
        inv_N = Arb(1 // N)
        (sqrt(λ_approx_div_λ(inv_N)) - 1) / inv_N^3
    end
    res1 = foldl(Arblib.union, values)

    # Compute enclosure for N >= N_max
    res2 = fx_div_x(Arb((0, 1 // N_max)), 3) do inv_N
        sqrt(λ_approx_div_λ(inv_N)) - 1
    end

    return Arblib.union(res1, res2)
end

"""
    sqrt_λ_approx_div_λ_remainder_N6(N₀::Int)

Compute an enclosure of
```
(sqrt(λ_approx / λ) - 1 - 2zeta(3) / N^3 - (6 - λ) * zeta(5) / N^5) * N^6
```
that is valid for all `N >= N₀`.
"""
function sqrt_λ_approx_div_λ_remainder_N6(N₀::Int)
    N_max = 1000

    # Compute enclosure for N₀ <= N <= N_max - 1
    λ = λ_disc()
    values = map(N₀:(N_max-1)) do N
        inv_N = Arb(1 // N)
        (
            sqrt(λ_approx_div_λ(inv_N)) - 1 - 2zeta(Arb(3)) * inv_N^3 -
            (6 - λ) * zeta(Arb(5)) * inv_N^5
        ) / inv_N^6
    end
    res1 = foldl(Arblib.union, values)

    # Compute enclosure for N >= N_max
    res2 = fx_div_x(Arb((0, 1 // N_max)), 6, force = true) do inv_N
        sqrt(λ_approx_div_λ(inv_N)) - 1 - 2zeta(Arb(3)) * inv_N^3 -
        (6 - λ) * zeta(Arb(5)) * inv_N^5
    end

    return Arblib.union(res1, res2)
end

"""
    _c_N(inv_N::Union{Arb,ArbSeries})

Compute an enclosure of `c_N`, given by
```
sqrt((gamma(1 - 1 / N)^2 * gamma(1 + 2 / N)) / (gamma(1 + 1 / N)^2 * gamma(1 - 2 / N)))
```
Note that this takes as input `inv(N)` and not `N`.

The computation is done using the formulation
```
rgamma(1 + 1 / N) / rgamma(1 - 1 / N) * sqrt(rgamma(1 - 2 / N) / rgamma(1 + 2 / N))
```
which is slightly more efficient. Here `rgamma(z) = 1 / gamma(z)` is
the reciprocal gamma function.
"""
_c_N(inv_N::Union{Arb,ArbSeries}) =
    rgamma(1 + inv_N) / rgamma(1 - inv_N) * sqrt(rgamma(1 - 2inv_N) / rgamma(1 + 2inv_N))

"""
    c_N(N₀::Int)

Compute an enclosure of `c_N` that is valid for all `N >= N₀`.
"""
c_N(N₀::Int) = 1 + c_N_remainder_mul_N3(N₀) * Arb((0, 1 // N₀))^3

"""
    c_N_remainder_mul_N3(N₀::Int)

Compute an enclosure of
```
(c_N - 1) * N^3
```
that is valid for all `N >= N₀`.
"""
function c_N_remainder_mul_N3(N₀::Int)
    N_max = 10000

    # Compute enclosure for N₀ <= N <= N_max - 1
    values = map(N₀:(N_max-1)) do N
        (_c_N(Arb(1 // N)) - 1) * N^3
    end
    res1 = foldl(Arblib.union, values)

    # Compute enclosure for N >= N_max
    res2 = fx_div_x(Arb((0, 1 // N_max)), 3, force = true) do inv_N
        _c_N(inv_N) - 1
    end

    return Arblib.union(res1, res2)
end

"""
    c_N_remainder_mul_N6(N₀::Int)

Compute an enclosure of
```
(c_N - 1 + 2zeta(3) / N^3 + 6zeta(5) / N^5) * N^6
```
that is valid for all `N >= N₀`.
"""
function c_N_remainder_mul_N6(N₀::Int)
    N_max = 20000

    # Compute enclosure for N₀ <= N <= N_max - 1
    values = map(N₀:(N_max-1)) do N
        inv_N = Arb(1 // N)
        (_c_N(inv_N) - 1 + 2zeta(Arb(3)) * inv_N^3 + 6zeta(Arb(5)) * inv_N^5) / inv_N^6
    end
    res1 = foldl(Arblib.union, values)

    # Compute enclosure for N >= N_max
    res2 = fx_div_x(Arb((0, 1 // N_max)), 6, force = true) do inv_N
        _c_N(inv_N) - 1 + 2zeta(Arb(3)) * inv_N^3 + 6zeta(Arb(5)) * inv_N^5
    end

    return Arblib.union(res1, res2)
end

function k_inv_N(inv_N::Union{Arb,ArbSeries})
    λ = λ_disc()
    R = Arb(R_inner)
    E_I =
        Arb(C_E_I_1) * inv_N +
        Arb(C_E_I_2) * inv_N^2 +
        Arb(C_E_I_3) * inv_N^3 +
        Arb(C_E_I_4) * inv_N^4
    a₀ = _c_N(inv_N) / (sqrt(λ) * besselj1(sqrt(λ)))

    return (
        2Arb(π) * a₀^2 * R^2 / 2 *
        (besselj0(R * sqrt(λ_approx(inv_N)))^2 + besselj1(R * sqrt(λ_approx(inv_N)))^2) -
        4Arb(π) * a₀ * E_I * R^2 / 2 *
        hypgeom0f1_regularized(Arb(2), -Arb(1 // 4) * R^2 * λ_approx(inv_N)) - R^2 * E_I^2
    )
end

function ϵ_prime(inv_N::Union{Arb,ArbSeries})
    ϵ = Arb(C_ε_6) * inv_N^6 + Arb(C_ε_7) * inv_N^7 + Arb(C_ε_8) * inv_N^8

    sqrt(Arb(π)) * ϵ / sqrt(k_inv_N(inv_N))
end
