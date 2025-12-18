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

Compute an enclosure `λ_app` from Equation REF(11) in the paper, given
by
```
λ * (1 + 4zeta(3) / N^3 + (12 - 2λ) * zeta(5) / N^5)
```
Note that this takes as input `inv(N)` and not `N`.
"""
λ_app(inv_N::Union{Arb,ArbSeries}) = λ_disc() * λ_app_div_λ(inv_N)

"""
    c_N(inv_N::Union{Arb,ArbSeries})

Compute an enclosure of `c_N` from Equation REF(6) in the paper, given
by
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
    F_N_model_remainder(N₀::Int, z::Acb)

Function for computing the remainder term used in [`F_N_model`](@ref).

It is an enclosure of the sixth derivative of `F_N(z)` in terms of
`inv(N)` on the interval ``[0, 1 / N₀]``, divided by `factorial(6)`.
The derivative is computed by moving it inside the integral expression
for `F_N(z)` and then enclosing the resulting integral.
"""
function F_N_model_remainder(N₀::Int, z::Acb)
    inv_N = Arb((0, 1 // N₀))

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

        # This enclosure is computed by explicitly computing the
        # derivative and writing it as a polynomial in log(t). The
        # integral is then split into terms, where all the factors
        # for the log(t) terms are bounded and can be factored out
        # from the integral.
        Arb((0, 1)) * (
            inv_N * powm1_div_t * integral_log(6, a) -
            6(2inv_N * (1 - t * z)^(-2inv_N) * log_1mtz_div_t - powm1_div_t) *
            integral_log(5, a) +
            60(1 - t * z)^(-2inv_N) *
            log_1mtz_div_t *
            (inv_N * log(1 - t * z) - 1) *
            integral_log(4, a) -
            80(2inv_N * log(1 - t * z) - 3) *
            log(1 - t * z) *
            (1 - t * z)^(-2inv_N) *
            log_1mtz_div_t *
            integral_log(3, a) +
            240(inv_N * log(1 - t * z) - 2) *
            log(1 - t * z)^2 *
            (1 - t * z)^(-2inv_N) *
            log_1mtz_div_t *
            integral_log(2, a) -
            96(2inv_N * log(1 - t * z) - 5) *
            log(1 - t * z)^3 *
            (1 - t * z)^(-2inv_N) *
            log_1mtz_div_t *
            integral_log(1, a) +
            64(inv_N * log(1 - t * z) - 3) *
            log(1 - t * z)^4 *
            (1 - t * z)^(-2inv_N) *
            log_1mtz_div_t *
            integral_log(0, a)
        )
    end

    # Integrate from a to b
    function integrand(t; analytic::Bool)
        if analytic
            # Check if t overlaps the branch cut. The branch cut
            # is for either t or 1 - tz lying on the negative real
            # axis. For the special case that n = 2 there is no
            # branch cut for t on the negative real axis.

            if Arblib.contains_nonpositive(Arblib.realref(t)) &&
               Arblib.contains_zero(Arblib.imagref(t))
                return indeterminate(t)
            end

            # Check if 1 - t * z overlaps the non-positive real axis.
            one_m_tz = 1 - t * z
            if Arblib.contains_nonpositive(Arblib.realref(one_m_tz)) &&
               Arblib.contains_zero(Arblib.imagref(one_m_tz))
                return indeterminate(t)
            end
        end

        ArbExtras.derivative_function(6) do inv_N
            inv_N * t^inv_N * ((1 - t * z)^-2inv_N - 1)
        end(inv_N) / t
    end

    remainder_a_b = Arblib.integrate(
        integrand,
        a,
        b,
        check_analytic = true,
        atol = radius(abs(remainder_0_a)) / 2,
        warn_on_no_convergence = false,
    )

    # Integrate from b to 1
    remainder_b_1 = if isone(b)
        zero(remainder_a_b) # Nothing to integrate
    else
        # Factor out bounds for log(t) as well as t^(inv_N - 1)
        # from the explicit integrands. The resulting integrals
        # can then be computed explicitly using
        # integral_logpow_1mtz.

        # Verify that the real and imaginary parts of log(1 - t *
        # z)^m * (1 - t * z)^-2inv_N don't change sign on the
        # interval, for 1 <= m <= 6.
        let t = Arb((b, 1))
            C = Arblib.abs_ubound(Arb, 1 - t * z)
            C < 1 || return indeterminate(z)

            # We are now ensured that log(1 - t * z) lies inside a
            # strip with real part (-∞, log(C)) and imaginary part (0,
            # π), (-π, 0) or [0, π] depending on weather imag(z) is
            # positive, negative or zero.

            # From the paper we then have that log(1 - t * z)^m *
            # (1 - t * z)^-2inv_N is contained in a single
            # quadrant if m * θ + 2inv_N * π <= π / 2. Where
            θ = atan(Arb(π), abs(log(C)))
            6θ + 2inv_N * π <= Arb(π) / 2 || return indeterminate(z)
        end

        let t = Arb((b, 1))
            t^(inv_N - 1) * (
                64inv_N * integral_logpow_1mtz(z, 6, -2inv_N, b) -
                192(inv_N * log(t) + 1) * integral_logpow_1mtz(z, 5, -2inv_N, b) +
                240(inv_N * log(t) + 2) * log(t) * integral_logpow_1mtz(z, 4, -2inv_N, b) -
                160(inv_N * log(t) + 3) *
                log(t)^2 *
                integral_logpow_1mtz(z, 3, -2inv_N, b) +
                60(inv_N * log(t) + 4) * log(t)^3 * integral_logpow_1mtz(z, 2, -2inv_N, b) -
                12(inv_N * log(t) + 5) * log(t)^4 * integral_logpow_1mtz(z, 1, -2inv_N, b) +
                (6 + inv_N * log(t)) *
                log(t)^5 *
                (integral_logpow_1mtz(z, 0, -2inv_N, b) - (1 - b))
            )
        end
    end

    (remainder_0_a + remainder_a_b + remainder_b_1) / factorial(6)
end

"""
    F_N_model(N₀::Int, z::Acb)

Compute an [`AcbTaylorModel`](@ref) of `F_N(z)` from Equation REF(7) in
the paper. The Taylor model is computed to degree 5 in terms of
`inv(N)` and is valid for all `inv(N)` in the interval ``[0,
inv(N₀)]``.

The details for the implementation are discussed in Appendix REF(B.4)
in the paper.
"""
function F_N_model(N₀::Int, z::Acb)
    return AcbTaylorModel(
        AcbSeries([1, 0, S(2, z), S(3, z), S(4, z), S(5, z), F_N_model_remainder(N₀, z)]),
        Arb((0, 1 // N₀)),
        Arb(0),
    )
end

"""
    epsilon(inv_N::Union{Arb,ArbSeries})

Compute `ε(N)` coming from Equation REF(23) in the paper. Note that
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
               inv_N * g_d2_1 * Arb(C_b_3) * Arb(C_T_4) +
               inv_N^2 * g_d2_1 / 2 * Arb(C_T_4)^2
           ) +
           inv_N^6 * (
        Arb(C_I_1_4) +
        Arb(C_I_2_3) +
        Arb(C_I_3_2) +
        Arb(C_I_K_1) +
        inv_N * (Arb(C_I_2_4) + Arb(C_I_3_3) + Arb(C_I_K_2)) +
        inv_N^2 * (Arb(C_I_3_4) + Arb(C_I_K_3)) +
        inv_N^3 * Arb(C_I_K_4)
    )
end

"""
    eta(inv_N::Union{Arb,ArbSeries})

Compute `η(N)` coming from Equation REF(26) in the paper. Note that
this takes as input `inv(N)` and not `N`.
"""
function eta(inv_N::Union{Arb,ArbSeries})
    λ = λ_disc()
    R = Arb(R_inner)
    E_I =
        Arb(C_V_1) * inv_N +
        Arb(C_V_2) * inv_N^2 +
        Arb(C_V_3) * inv_N^3 +
        Arb(C_V_4) * inv_N^4
    a₀ = c_N(inv_N) / (sqrt(λ) * besselj1(sqrt(λ)))

    return (
        Arb(π) *
        a₀^2 *
        R^2 *
        (besselj0(R * sqrt(λ_app(inv_N)))^2 + besselj1(R * sqrt(λ_app(inv_N)))^2) -
        4Arb(π) * a₀ * E_I * R / sqrt(λ_app(inv_N)) * besselj1(R * sqrt(λ_app(inv_N))) -
        π * R^2 * E_I^2
    )
end

"""
    epsilon_hat(inv_N::Union{Arb,ArbSeries})

Compute `hat{ε}(N)` coming from Equation REF(27) in the paper. Note
that this takes as input `inv(N)` and not `N`.
"""
function epsilon_hat(inv_N::Union{Arb,ArbSeries})
    return sqrt(Arb(π)) * epsilon(inv_N) / sqrt(eta(inv_N))
end
