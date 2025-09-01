# This file contains the implementation the functions that are bounded
# in Lemma 2.6.

function F_N(inv_N::Arb, z::Acb)
    return 1 + inv_N * F_N_sub_1_mul_N(inv_N, z)
end

function F_N_sub_1_mul_N(inv_N::Arb, z::Acb)
    a = Arb(1e-8)
    b = Arblib.contains(z, Acb(1)) ? Arb(1) - 1e-8 : Arb(1)

    # Integrate from 0 to a
    res_0_a = a * Arb((0, 1)) * fx_div_x(Acb(Arblib.union(Arb(0), a))) do t
        (1 - t * z)^-2inv_N - 1
    end

    # Integrate from a to b
    res_a_b = Arblib.integrate(a, b, warn_on_no_convergence = false) do t
        t^inv_N * ((1 - t * z)^-2inv_N - 1) / t
    end

    # Integrate from b to 1
    res_b_1 = if isone(b)
        zero(res_a_b) # Nothing to integrate
    else
        let t = Arblib.union(b, Arb(1))
            t^inv_N / t * (
                ((1 - b * z)^(1 - 2inv_N) - pow(1 - z, 1 - 2inv_N)) / (z * (1 - 2inv_N)) - (1 - b)
            )
        end
    end

    return res_0_a + res_a_b + res_b_1
end

"""
    F_N_model(N₀::Int, z::Acb)

Compute a degree 6 [`ArbTaylorModel`](@ref) of [`F_N`](@ref) in
`inv(N)` that is valid for all `N >= N₀`.
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
        AcbSeries([
            1,
            0,
            S_integral(2, z),
            S_integral(3, z),
            S_integral(4, z),
            S_integral(5, z),
            F_N_remainder,
        ]),
        inv_N,
        Arb(0),
    )
end

function b_2(z::Acb)
    return real(S_integral(2, z))
end

function b_3(z::Acb)
    return real(S_integral(3, z))
end

function b_4(z::Acb)
    S2_imag = imag(S_integral(2, z))
    S4_real = real(S_integral(4, z))

    return S2_imag^2 / 2 + S4_real
end

function b_5(z::Acb)
    S2_imag = imag(S_integral(2, z))
    S3_imag = imag(S_integral(3, z))
    S5_real = real(S_integral(5, z))

    return S2_imag * S3_imag + S5_real - λ_disc() * zeta(Arb(5))
end

function T_6(N₀::Int)
    inv_N = Arb((0, 1 // N₀))

    # It is really important to get a good enclosure of ρ, so we
    # manual bisect to get better enclosures
    M1 = truncate(
        ArbTaylorModel(inv_N, Arb(0), degree = 10) do inv_N
            if iswide(inv_N[0])
                # Bisect in inv_N[0]. Note that all other coefficients
                # are always exact.
                parts = map(
                    ArbExtras.bisect_interval_recursive(
                        Arblib.getinterval(inv_N[0])...,
                        10,
                    ),
                ) do inv_N0_part
                    inv_N_part = copy(inv_N)
                    inv_N_part[0] = inv_N0_part
                    sqrt(λ_approx_div_λ(inv_N_part)) * _c_N(inv_N_part)
                end
                foldl(Arblib.union, parts)
            else
                sqrt(λ_approx_div_λ(inv_N)) * _c_N(inv_N)
            end
        end,
        degree = 5,
    )

    return z::Acb -> begin
        M2 = abs(F_N_model(N₀, z))

        return (M1 * M2).p[6]
    end
end
