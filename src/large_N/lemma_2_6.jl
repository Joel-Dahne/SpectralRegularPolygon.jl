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

        # Integrate from a to b
        remainder_a_b = Arblib.integrate(a, b) do t
            ArbExtras.derivative_function(6) do inv_N
                inv_N * t^inv_N * ((1 - t * z)^-2inv_N - 1) / t
            end(inv_N)
        end

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

function abs_F_N_remainder_N6(N₀::Int, z::Acb)
    inv_N = Arb((0, 1 // N₀))

    # Step 1: Compute a Taylor model of F_N(z) in N^-1 with remainder
    # term of degree 6.

    # Step 1.1: Compute bound on remainder term. This is an enclosure
    # of the sixth derivative of F_N(z) in N, divided by factorial(6).
    # To compute the derivative we move the derivative inside the
    # integral.
    F_N_remainder = let
        a = Arb(1e-8)
        b = Arblib.contains(z, Acb(1)) ? Arb(1) - 1e-8 : Arb(1)

        # Integrate from a to b
        remainder_a_b = Arblib.integrate(a, b) do t
            ArbExtras.derivative_function(6) do inv_N
                inv_N * t^inv_N * ((1 - t * z)^-2inv_N - 1) / t
            end(inv_N)
        end

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

    # Step 1.2: Compute Taylor series
    F_N_model = AcbTaylorModel(
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
    # Step 2: Compute a Taylor model of abs(F_N(z))
    abs_F_N_model = abs(F_N_model)

    return real(abs_F_N_model.p[end]) # Return bound on remainder
end

# b_2(z) = S(2, z) + S(2, conj(z))) / 2
# NOTE: This assumes that abs(z) == 1
function b_2(z::Acb)
    return real(S_integral(2, z))

    if Arblib.contains(z, Acb(1, 0))
        # Input overlaps branch cut
        return real(S_unitdisc(2, z))
    else
        return real(mean_value_theorem_bound(z -> S(2, z), z))
    end
end

# b_3(z) = (S(3, z) + S(3, conj(z))) / 2
function b_3(z::Acb)
    return real(S_integral(3, z))
    return real(mean_value_theorem_bound(z -> S(3, z), z))
end

# b_4(z) = -(S(2, z) - S(2, conj(z)))^2 / 8 + (S(4, z) + S(4, conj(z))) / 2
function b_4(z::Acb)
    S2_imag = imag(S_integral(2, z))
    S4_real = real(S_integral(4, z))

    #S2_imag = imag(mean_value_theorem_bound(z -> S(2, z), z))
    #S4_real = real(mean_value_theorem_bound(z -> S(4, z), z))

    return S2_imag^2 / 2 + S4_real
end

#b_5(z) =
#    (
#        -(S(2, z) - S(2, conj(z))) * (S(3, z) - S(3, conj(z))) +
#        2(S(5, z) + S(5, conj(z)) - 2λ_disc() * zeta(Arb(5)))
#    ) / 4
function b_5(z::Acb)
    S2_imag = imag(S_integral(2, z))
    S3_imag = imag(S_integral(3, z))
    S5_real = real(S_integral(5, z))

    #S2_imag = imag(mean_value_theorem_bound(z -> S(2, z), z))
    #S3_imag = imag(mean_value_theorem_bound(z -> S(3, z), z))
    #S5_real = real(mean_value_theorem_bound(z -> S(5, z), z))

    return S2_imag * S3_imag + S5_real - λ_disc() * zeta(Arb(5))
end
_b_5_real(z) = S(5, z) - λ_disc() * zeta(Arb(5))
_b_5_imag(z) = imag(S(2, z)) * imag(S(3, z))

function T_6_bound(N₀::Int)
    # Compute all parts not depending on z
    inv_N = Arb((0, 1 // N₀))

    # Compute factors with removable singularities

    # (sqrt(λ_approx / λ) - 1 - 2zeta(3) / N^3 - (6 - λ) * zeta(5) / N^5) * N^6
    removable_1 = sqrt_λ_approx_div_λ_remainder_N6(N₀)

    # (c_N - 1 + 2zeta(3) / N^3 + 6zeta(5) / N^5) * N^6
    removable_2 = c_N_remainder_mul_N6(N₀)

    # (sqrt(λ_approx / λ) - 1) * N^3
    removable_3 = sqrt_λ_approx_div_λ_remainder_N3(N₀)

    # Compute factors not depending on z

    # Factor for abs(F_N(z) - 1) * N
    factor_term2 = abs(λ_disc() * zeta(Arb(5)))

    # Factor for abs(F_N(z))
    factor_term3 = abs(
        removable_1 +
        sqrt(λ_approx_div_λ(inv_N)) * removable_2 +
        (-2zeta(Arb(3)) + 6zeta(Arb(5)) * inv_N^2) * removable_3,
    )

    return z -> let
        # Enclosure of
        # (abs(F_N(z) - 1 - b_2(z) / N^2 - b_3(z) / N^3 - b_4(z) / N^4 - XXX)) * N^6
        term1 = abs_F_N_remainder_N6(N₀, z)

        # Enclosure of (F_N(z) - 1) * N
        F_N_sub_1_mul_N_enclosure = F_N_sub_1_mul_N(inv_N, z)

        # Bound of (abs(F_N(z)) - 1) * N
        term2_part = abs(F_N_sub_1_mul_N(inv_N, z))

        term2 = term2_part * factor_term2

        # Enclosure of abs(F_N(z))
        term3_part = abs(1 + inv_N * F_N_sub_1_mul_N_enclosure)

        term3 = term3_part * factor_term3

        return term1 + term2 + term3
    end
end
