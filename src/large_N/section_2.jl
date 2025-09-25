"""
    λ_disc()

Compute an enclosure of the first eigenvalue of the unit disc.

It is computed as the square of the first positive root of `besselj0`.
"""
function λ_disc()
    # Isolate all roots on the interval [0, 2.5]
    roots, flags = ArbExtras.isolate_roots(besselj0, Arf(0), Arf(2.5))

    # Verify that there is exactly one root, this proves that it is
    # the first one.
    length(flags) == 1 && flags[1] || throw(ErrorException("could not isolate root"))

    # Refine the enclosure
    return ArbExtras.refine_root(besselj0, Arb(only(roots)))^2
end

"""
    λ_app_div_λ(inv_N::Union{Arb,ArbSeries})

Compute an enclosure of [`λ_app`](@ref) divided by `λ`, given by
```
1 + 4zeta(3) / N^3 + (12 - 2λ) * zeta(5) / N^5
```
Note that this takes as input `inv(N)` and not `N`.
"""
λ_app_div_λ(inv_N::Union{Arb,ArbSeries}) =
    1 + 4zeta(Arb(3)) * inv_N^3 + (12 - 2λ_disc()) * zeta(Arb(5)) * inv_N^5

"""
    λ_app(inv_N::Union{Arb,ArbSeries})

Compute an enclosure of
```
λ * (1 + 4zeta(3) / N^3 + (12 - 2λ) * zeta(5) / N^5)
```
Note that this takes as input `inv(N)` and not `N`.
"""
λ_app(inv_N::Union{Arb,ArbSeries}) = λ_disc() * λ_app_div_λ(inv_N)

"""
    c_N(inv_N::Union{Arb,ArbSeries})

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
function c_N(inv_N::Union{Arb,ArbSeries})
    if inv_N isa ArbSeries && iswide(inv_N[0])
        # T_6 and K_4 both require accurate enclosures of high degree
        # expansions of c_N. Direct evaluation gives fairly poor
        # bounds when inv_N is moderately wide. For that reason we
        # manually bisect the constant part of inv_N and take the
        # union of the results. Note that this only makes sense
        # because in practice all coefficients but the first one of
        # inv_N are exact.

        # Constant part of the expansion
        inv_N_0 = inv_N[0]
        # Bisect the constant part recursively 10 times
        inv_N_0_parts = ArbExtras.bisect_interval_recursive(getinterval(inv_N_0)..., 10)


        c_N_parts = map(inv_N_0_parts) do inv_N_0_part
            # Set the expansion
            inv_N_part = copy(inv_N)
            inv_N_part[0] = inv_N_0_part
            # Evaluate on current part
            rgamma(1 + inv_N_part) / rgamma(1 - inv_N_part) *
            sqrt(rgamma(1 - 2inv_N_part) / rgamma(1 + 2inv_N_part))
        end

        # Return union of all parts
        return foldl(Arblib.union, c_N_parts)
    else
        return rgamma(1 + inv_N) / rgamma(1 - inv_N) *
               sqrt(rgamma(1 - 2inv_N) / rgamma(1 + 2inv_N))
    end
end

"""
    F_N_model(N₀::Int, z::Acb)

Compute an [`ArbTaylorModel`](@ref) of [`F_N`](@ref) in `inv(N)` that
is valid for all `N >= N₀`. It is computed with a remainder term of
degree 6.

TODO: Write documentation for this.
"""
function F_N_model(N₀::Int, z::Acb)
    inv_N = Arb((0, 1 // N₀))

    # Compute bound on remainder term. This is an enclosure of the
    # sixth derivative of F_N(z) in N, divided by factorial(6). To
    # compute the derivative we move the derivative inside the
    # integral.
    F_N_remainder = let
        a = Arb(1e-8)
        b = Arblib.contains(z, Acb(1)) ? Arb(1) - 1e-8 : Arb(1)

        # Integrate from 0 to a
        remainder_0_a = let t = Arb((0, a))
            # Enclosure of log(1 - t * z) / t
            log_1mtz_div_t = fx_div_x(Acb(t)) do t
                log(1 - t * z)
            end
            # Enclosure of ((1 - t * z)^(-2inv_N) - 1) / t
            powm1_div_t = fx_div_x(Acb(t)) do t
                (1 - t * z)^(-2inv_N) - 1
            end

            # PROVE: Write documentation for this. It is based on
            # explicitly computing the 6th derivative and then
            # factoring out all bounded terms from the integral. This
            # leaves only integrals of power of logarithms, which are
            # computed explicitly using integral_logm.
            Arb((0, 1)) * (
                integral_log(6, a) * inv_N * powm1_div_t +
                integral_log(5, a) *
                (6powm1_div_t - 12inv_N * (1 - t * z)^(-2inv_N) * log_1mtz_div_t) +
                integral_log(4, a) *
                60(1 - t * z)^(-2inv_N) *
                log_1mtz_div_t *
                (-1 + inv_N * log(1 - t * z)) +
                integral_log(3, a) *
                80(1 - t * z)^(-2inv_N) *
                log_1mtz_div_t *
                log(1 - t * z) *
                (3 - 2inv_N * log(1 - t * z)) +
                integral_log(2, a) *
                240(1 - t * z)^(-2inv_N) *
                log_1mtz_div_t *
                log(1 - t * z)^2 *
                (-2 + inv_N * log(1 - t * z)) +
                integral_log(1, a) *
                96(1 - t * z)^(-2inv_N) *
                log_1mtz_div_t *
                log(1 - t * z)^3 *
                (5 - 2inv_N * log(1 - t * z)) +
                integral_log(0, a) *
                64(1 - t * z)^(-2inv_N) *
                log_1mtz_div_t *
                log(1 - t * z)^4 *
                (-3 + inv_N * log(1 - t * z))
            )
        end

        # Integrate from a to b
        remainder_a_b = Arblib.integrate(
            a,
            b,
            atol = Arblib.radius(abs(remainder_0_a)) / 2,
            warn_on_no_convergence = false,
        ) do t
            ArbExtras.derivative_function(6) do inv_N
                inv_N * t^inv_N * ((1 - t * z)^-2inv_N - 1)
            end(inv_N) / t
        end

        # Integrate from b to 1
        remainder_b_1 = if isone(b)
            zero(remainder_a_b) # Nothing to integrate
        else
            let t = Arb((b, 1))
                # PROVE: Write documentation for this. It is based on
                # explicitly computing the 6th derivative and then
                # factoring out all bounded terms from the integral.
                # This leaves only integrals of the form log(1 - t *
                # z)^m * (1 - t * z)^(-2inv_N), which are explicitly
                # computed using integral_logpow_1mtz.
                t^(-1 + inv_N) * (
                    log(t)^5 *
                    (6 + inv_N * log(t)) *
                    (integral_logpow_1mtz(z, 0, -2inv_N, b) - (1 - b)) -
                    12log(t)^4 *
                    (5 + inv_N * log(t)) *
                    integral_logpow_1mtz(z, 1, -2inv_N, b) +
                    60log(t)^3 *
                    (4 + inv_N * log(t)) *
                    integral_logpow_1mtz(z, 2, -2inv_N, b) -
                    160log(t)^2 *
                    (3 + inv_N * log(t)) *
                    integral_logpow_1mtz(z, 3, -2inv_N, b) +
                    240log(t) *
                    (2 + inv_N * log(t)) *
                    integral_logpow_1mtz(z, 4, -2inv_N, b) -
                    192(1 + inv_N * log(t)) * integral_logpow_1mtz(z, 5, -2inv_N, b) +
                    64inv_N * integral_logpow_1mtz(z, 6, -2inv_N, b)
                )
            end
        end

        (remainder_0_a + remainder_a_b + remainder_b_1) / factorial(6)
    end

    return AcbTaylorModel(
        AcbSeries([1, 0, S(2, z), S(3, z), S(4, z), S(5, z), F_N_remainder]),
        inv_N,
        Arb(0),
    )
end

"""
    epsilon(inv_N::Union{Arb,ArbSeries})

Compute `ε(N)` coming from Equation REF(16) in the paper. Note that
this takes as input `inv(N)` and not `N`.
"""
function epsilon(inv_N::Union{Arb,ArbSeries})
    # Bound of g'(1), which is the same as that for g''(1)
    λ = λ_disc()
    g_d_1 = abs(sqrt(λ) * besselj1(sqrt(λ)))
    g_d2_1 = g_d_1

    return Arb(C_a_0) *
           inv_N^6 *
           (
               g_d_1 * Arb(C_T_6) +
               g_d2_1 / 2 * Arb(C_b_3)^2 +
               g_d2_1 * Arb(C_b_2) * Arb(C_T_4) +
               Arb(C_T_2)^3 * Arb(C_gd3) / 6 +
               2inv_N * Arb(C_b_3) * Arb(C_T_4) +
               inv_N^2 * Arb(C_T_4)^2
           ) + inv_N^6 * (Arb(C_I_1_4) + Arb(C_I_2_3) + Arb(C_I_3_2) + Arb(C_I_4_1))
end

"""
    k_inv_N(inv_N::Union{Arb,ArbSeries})

Compute `k(N)` coming from Equation REF(18) in the paper. Note that
this takes as input `inv(N)` and not `N`.
"""
function k_inv_N(inv_N::Union{Arb,ArbSeries})
    λ = λ_disc()
    R = Arb(R_inner)
    E_I =
        Arb(C_V_1) * inv_N +
        Arb(C_V_2) * inv_N^2 +
        Arb(C_V_3) * inv_N^3 +
        Arb(C_V_4) * inv_N^4
    a₀ = c_N(inv_N) / (sqrt(λ) * besselj1(sqrt(λ)))

    return (
        2Arb(π) * a₀^2 * R^2 / 2 *
        (besselj0(R * sqrt(λ_app(inv_N)))^2 + besselj1(R * sqrt(λ_app(inv_N)))^2) -
        4Arb(π) * a₀ * E_I * R^2 / 2 *
        hypgeom0f1_regularized(Arb(2), -Arb(1 // 4) * R^2 * λ_app(inv_N)) - R^2 * E_I^2
    )
end

"""
    epsilon_hat(inv_N::Union{Arb,ArbSeries})

Compute `hat{ε}` coming from Equation REF(19) in the paper. Note
that this takes as input `inv(N)` and not `N`.
"""
function epsilon_hat(inv_N::Union{Arb,ArbSeries})
    return sqrt(Arb(π)) * epsilon(inv_N) / sqrt(k_inv_N(inv_N))
end
