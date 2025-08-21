# PROVE: We need to verify that these are indeed the first zeros. For
# the first one we could do this with ArbExtras.isolate_roots. For the
# second one we have to explicitly handle the neighbourhood of zero,
# but that is fine. For the computational cost it would however be
# nice to not have to do this every time.
λ_disc() = ArbExtras.refine_root(besselj0, Arb((sqrt(Arf(5.7)), sqrt(Arf(5.8)))))^2
λ₂_disc() = ArbExtras.refine_root(besselj1, Arb((Arf(3.75), Arf(3.875))))^2

function λ_approx(inv_N::Union{Arb,ArbSeries})
    λ = λ_disc()
    return λ * (1 + 4zeta(Arb(3)) * inv_N^3 + (12 - 2λ) * zeta(Arb(5)) * inv_N^5)
end
λ_approx_div_λ(inv_N::Union{Arb,ArbSeries}) =
    1 + 4zeta(Arb(3)) * inv_N^3 + (12 - 2λ_disc()) * zeta(Arb(5)) * inv_N^5

# IMPROVE: We could cache this value
# Enclosure of (sqrt(λ_approx / λ) - 1) * N^3 for N >= N₀.
function sqrt_λ_approx_div_λ_remainder_N3(N₀::Int)
    N_max = 100
    λ = λ_disc()
    # Compute enclosure for N from N₀ to N_max
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

# IMPROVE: We could cache this value
# Enclosure of (sqrt(λ_approx / λ) - ...) * N^6 for N >= N₀,
# where ... denotes the first three terms in the expansion of
# sqrt(λ_approx / λ) in inv(N).
function sqrt_λ_approx_div_λ_remainder_N6(N₀::Int)
    N_max = 1000
    λ = λ_disc()
    # Compute enclosure for N from N₀ to N_max
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

# c_N but computed using rgamma instead of gamma since that gives
# better enclosures.
_c_N(inv_N::Union{Arb,ArbSeries}) = sqrt(
    (rgamma(1 + inv_N)^2 * rgamma(1 - 2inv_N)) / (rgamma(1 - inv_N)^2 * rgamma(1 + 2inv_N)),
)

# IMPROVE: We could cache this value
# Enclosure of (c_N - 1) * N^3 for N >= N₀
function c_N_remainder_mul_N3(N₀::Int)
    N_max = 10000
    # Compute enclosure for N from N₀ to N_max
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

# IMPROVE: We could cache this value
# Enclosure of (c_N - ...) * N^6 for N >= N₀, where ... denotes the
# first three terms in the expansion of c_N in inv(N).
function c_N_remainder_mul_N6(N₀::Int)
    N_max = 20000
    # Compute enclosure for N from N₀ to N_max
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

# Enclosure of C_N for N >= N₀
c_N(N₀::Int) = 1 + c_N_remainder_mul_N3(N₀) * Arb((0, 1 // N₀))^3

function F_N_sub_1_mul_N(inv_N::Arb, z::Acb)
    a = Arb(1e-8)
    b = Arblib.contains(z, Acb(1)) ? Arb(1) - 1e-8 : Arb(1)

    # Integrate from a to b
    res_a_b = Arblib.integrate(a, b) do t
        t^inv_N * ((1 - t * z)^-2inv_N - 1) / t
    end

    # Integrate from 0 to a
    res_0_a = a * Arb((0, 1)) * fx_div_x(Acb(Arblib.union(Arb(0), a))) do t
        (1 - t * z)^-2inv_N - 1
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

function F_N(inv_N::Arb, z::Acb)
    return 1 + inv_N * F_N_sub_1_mul_N(inv_N, z)
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
    F_N_model = TaylorModel(
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
        #@show (term1, term2, term2)
        return term1 + term2 + term3
    end
end

function g_dw3(w)
    λ = λ_disc()
    return λ^(3 // 2) * (3besselj1(w * sqrt(λ) - besselj(oftype(λ, 3), w * sqrt(λ)))) / 4
end

V_1(z) = 2polylog(1, z)

function V_1_div_z_bound(zᵤ::Arb)
    0 < zᵤ < 1 || return indeterminate(zᵤ)
    return 2polylog_r_div_z_bound(1, zᵤ)
end

function V_2(z)
    λ = λ_disc()
    return (λ / 2 - 2) * polylog(2, z) + 4polylog_1_1(z)
end

function V_2_div_z_bound(zᵤ::Arb)
    0 < zᵤ < 1 || return indeterminate(zᵤ)
    λ = λ_disc()
    return abs(λ / 2 - 2) * polylog_r_div_z_bound(1, zᵤ) + 4polylog_r_div_z_bound(2, zᵤ)
end

function V_3(z)
    λ = λ_disc()
    return (λ^2 / 16 - λ + 2) * polylog(3, z) +
           (3λ - 12) * polylog_1_2(z) +
           (λ - 4) * polylog_2_1(z) - 8polylog_1_1_1(z)
end

function V_3_div_z_bound(zᵤ::Arb)
    0 < zᵤ < 1 || return indeterminate(zᵤ)
    λ = λ_disc()
    return abs(λ^2 / 16 - λ + 2) * polylog_r_div_z_bound(1, zᵤ) +
           abs(3λ - 12) * polylog_r_div_z_bound(2, zᵤ) +
           abs(λ - 4) * polylog_r_div_z_bound(2, zᵤ) +
           8polylog_r_div_z_bound(3, zᵤ)
end

function V_4(z)
    λ = λ_disc()
    return (λ^3 / 192 - λ^2 / 8 - λ / 2 - 2) * polylog(4, z) +
           (λ^2 / 8 - 2λ + 4) * polylog_3_1(z) +
           (λ^2 / 4 - 4λ + 12) * polylog_2_2(z) +
           (5λ^2 / 8 - 8λ + 28) * polylog_1_3(z) +
           (2λ - 8) * polylog_2_1_1(z) +
           (6λ - 24) * polylog_1_2_1(z) +
           (14λ - 56) * polylog_1_1_2(z) +
           16polylog_1_1_1_1(z) +
           2λ * zeta(Arb(3)) * polylog(1, z)
end

function V_4_div_z_bound(zᵤ::Arb)
    0 < zᵤ < 1 || return indeterminate(zᵤ)
    λ = λ_disc()
    return abs(λ^3 / 192 - λ^2 / 8 - λ / 2 - 2) * polylog_r_div_z_bound(1, zᵤ) +
           abs(λ^2 / 8 - 2λ + 4) * polylog_r_div_z_bound(2, zᵤ) +
           abs(λ^2 / 4 - 4λ + 12) * polylog_r_div_z_bound(2, zᵤ) +
           abs(5λ^2 / 8 - 8λ + 28) * polylog_r_div_z_bound(2, zᵤ) +
           abs(2λ - 8) * polylog_r_div_z_bound(3, zᵤ) +
           abs(6λ - 24) * polylog_r_div_z_bound(3, zᵤ) +
           abs(14λ - 56) * polylog_r_div_z_bound(3, zᵤ) +
           16polylog_r_div_z_bound(4, zᵤ) +
           abs(2λ * zeta(Arb(3))) * polylog_r_div_z_bound(1, zᵤ)
end

V(l::Int, z) =
    if l == 1
        V_1(z)
    elseif l == 2
        V_2(z)
    elseif l == 3
        V_3(z)
    elseif l == 4
        V_4(z)
    end

V_div_z_bound(l::Int, zᵤ::Arb) =
    if l == 1
        V_1_div_z_bound(zᵤ)
    elseif l == 2
        V_2_div_z_bound(zᵤ)
    elseif l == 3
        V_3_div_z_bound(zᵤ)
    elseif l == 4
        V_4_div_z_bound(zᵤ)
    end

# Return C s.t. abs(V_2(z)) <= C * abs(z)
function V_2_bound(z::Acb)
    (λ_disc() / 2 - 2) * polylog(2, z) + 2log(1 - z)^2
end

d_0(z, t) = one(z)
d_0_z_part(z) = one(z)
d_0_zt_part(z, t) = zero(z)

d_1(z, t) = λ_disc() / 4 * log(t / z)
d_1_z_part(z) = zero(z)
d_1_zt_part(z, t) = λ_disc() / 4 * log(t / z)

function d_2(z, t)
    λ = λ_disc()
    return λ^2 / 64 * log(t / z)^2 +
           λ / 8 * (log(t / z)^2 + 2S_integral(2, t) - 2S_integral(2, z))
end
function d_2_z_part(z)
    λ = λ_disc()
    return -λ / 4 * S_integral(2, z)
end
function d_2_zt_part(z, t)
    λ = λ_disc()
    return λ^2 / 64 * log(t / z)^2 + λ / 8 * (log(t / z)^2 + 2S_integral(2, t))
end

function d_3(z, t)
    λ = λ_disc()
    # NOTE: This is an alternative formulation that is SLIGHTLY better
    # than the version in the paper.
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + λ / 16 + 1 // 6) * log(t / z)^2 + (λ / 8 + 1) * S_integral(2, t) -
            λ / 8 * S_integral(2, z) + S_integral(2, conj(z))
        ) + S_integral(3, t) - S_integral(3, z)
    )
    # This is the version in the paper
    #return λ^3 / 2304 * log(t / z)^3 +
    #       λ^2 / 64 * log(t / z) * (log(t / z)^2 + 2S_integral(2, t) - 2S_integral(2, z)) +
    #       λ / 24 * (
    #           log(t / z)^3 +
    #           6log(t / z) * S_integral(2, t) +
    #           6log(t / z) * S_integral(2, conj(z)) +
    #           6S_integral(3, t) - 6S_integral(3, z)
    #       )
end
function d_3_z_part(z)
    λ = λ_disc()
    return -λ / 4 * S_integral(3, z)
end
function d_3_zt_part(z, t)
    λ = λ_disc()
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + λ / 16 + 1 // 6) * log(t / z)^2 + (λ / 8 + 1) * S_integral(2, t) -
            λ / 8 * S_integral(2, z) + S_integral(2, conj(z))
        ) + S_integral(3, t)
    )
end

# Part of d_3 analytic in z
function d_3_analytic_z(z, t)
    λ = λ_disc()
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + (1 // 16) * λ + 1 // 6) * log(t / z)^2 +
            ((1 // 8) * λ + 1) * S_integral(2, t) - (1 // 8) * λ * S_integral(2, z)
        ) + S_integral(3, t) - S_integral(3, z)
    )
end

# Part of d_3 analytic in conj(z)
function d_3_analytic_conj_z(z, t)
    λ = λ_disc()
    return λ / 4 * log(t / z) * S_integral(2, conj(z))
end

function integral_d_0(z::Acb, a::Arb)
    return a * z
end

function integral_d_1(z::Acb, a::Arb)
    λ = λ_disc()
    return λ / 4 * integral_log_z(1, z, a)
end

"""
    integral_d_2(z::Acb, a::Arb)

We write `d_2` as
```
(λ^2 / 64 + λ / 8) * log(t / z)^2 +
    λ / 4 * (S(2, t) - S(2, z))
```
The `S` factors with `t` we enclose on the entire interval.
We then integrate the log-terms with [`integral_log_z`](@ref).
"""
function integral_d_2(z::Acb, a::Arb)
    λ = λ_disc()
    t = Arblib.union(zero(z), a * z)
    return (λ^2 / 64 + λ / 8) * integral_log_z(2, z, a) +
           λ / 4 * (S_integral(2, t) - S_integral(2, z)) * a * z
end

"""
    integral_d_3(z::Acb, a::Arb)

We write `d_3` as
```
(λ^3 / 2304 + λ^2 / 64 + λ / 24) * log(t / z)^3 +
    (
        (λ^2 / 32 + λ / 4) * S(2, t) -
        (λ^2 / 32 * S(2, z) - λ / 4 * S(2, conj(z)))
    ) * log(t / z)
    λ / 4 * (S(3, t) - S(3, z))
```
The `S` factors with `t` we enclose on the entire interval.
We then integrate the log-terms with [`integral_log_z`](@ref).
"""
function integral_d_3(z::Acb, a::Arb)
    λ = λ_disc()
    t = Arblib.union(zero(z), a * z)
    return (λ^3 / 2304 + λ^2 / 64 + λ / 24) * integral_log_z(3, z, a) +
           (
               (λ^2 / 32 + λ / 4) * S_integral(2, t) -
               (λ^2 / 32 * S_integral(2, z) - λ / 4 * S_integral(2, conj(z)))
           ) * integral_log_z(1, z, a) +
           λ / 4 * (S_integral(3, t) - S_integral(3, z)) * a * z
end

d(k::Int, z, t) =
    if k == 0
        d_0(z, t)
    elseif k == 1
        d_1(z, t)
    elseif k == 2
        d_2(z, t)
    elseif k == 3
        d_3(z, t)
    end

# Part of d only depending on z
d_z_part(k::Int, z) =
    if k == 0
        d_0_z_part(z)
    elseif k == 1
        d_1_z_part(z)
    elseif k == 2
        d_2_z_part(z)
    elseif k == 3
        d_3_z_part(z)
    end

# Part of d only depending on z and t
d_zt_part(k::Int, z, t) =
    if k == 0
        d_0_zt_part(z, t)
    elseif k == 1
        d_1_zt_part(z, t)
    elseif k == 2
        d_2_zt_part(z, t)
    elseif k == 3
        d_3_zt_part(z, t)
    end

"""
    integral_d(k::Int, z::Acb, a::Arb)

Compute the integral of `d(l, z, t)` from `0` to `a * z`.
"""
integral_d(k::Int, z::Acb, a::Arb) =
    if k == 0
        integral_d_0(z, a)
    elseif k == 1
        integral_d_1(z, a)
    elseif k == 2
        integral_d_2(z, a)
    elseif k == 3
        integral_d_3(z, a)
    end

"""
    integral_d_0_V_1(z::Acb)

Integrating
```
1 / t * d(0, z, t) * V(1, t) = -2log(1 - t) / t
```
Gives us
```
2polylog(2, t)
```
Which from `0` to `z` gives us
```
2polylog(2, z)
```
"""
function integral_d_0_V_1(z::Acb)
    return 2polylog2(z)
end

"""
    integral_d_1_V_1(z::Acb)

Integrating
```
1 / t * d(1, z, t) * V(1, t) = -λ / 2 * log(t / z) * log(1 - t) / t
```
Gives us
```
λ / 2 * (log(t / z) * polylog(2, t) + polylog(3, t))
```
Which from `0` to `z` gives us
```
λ / 2 * polylog(3, z)
```
"""
function integral_d_1_V_1(z::Acb)
    return λ_disc() / 2 * polylog(3, z)
end

function integral_d_k_V_l(k::Int, l::Int, z::Acb)
    if k == 0 && l == 1
        return integral_d_0_V_1(z)
    elseif k == 1 && l == 1
        return integral_d_1_V_1(z)
    end

    # NEXT: Check the enclosures don't convert for (k, l) = (3, 2)

    a = Arb(1e-4)
    az = a * z

    # Integrate from 0 to a * z
    res_0_az = let
        # Factor out enclosure of V(l, t) / t and integrate d(k, z, t)

        # Compute enclosure of V(l, t) / t for t in [0, a * z]
        # IMPROVE: Improving this enclosure would allow us to take a
        # larger a, which should make the integration much faster.
        azᵤ = abs_ubound(Arb, az)
        V_div_t = add_error(Acb(0), V_div_z_bound(l, azᵤ))

        (V_div_t) * (integral_d(k, z, a))
    end

    # Integrate from a * z to z
    # TODO: We need to verify analyticity for this to be correct,
    # which we might not have.
    res_az_z = let
        # Arblib.integrate doesn't handle wide integration limits very
        # well. For that reason we integrate from the midpoint of the
        # endpoints if z is wide, and add the remaining part later.
        az_thin, z_thin = if iswide(z)
            midpoint(Acb, az), midpoint(Acb, z)
        else
            az, z
        end

        # Part of d(k, z, t) only depending on z, so we can factor it
        # out.

        res_az_thin_z_thin_part_1 =
            d_z_part(k, z) * Arblib.integrate(
                az_thin,
                z_thin,
                atol = 1e-6,
                opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
                warn_on_no_convergence = false,
            ) do t
                V(l, t) / t
            end

        # Part of d(k, z, t) depending on both z and t.
        res_az_thin_z_thin_part_2 =
            Arblib.integrate(
                az_thin,
                z_thin,
                atol = 1e-6,
                opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
                warn_on_no_convergence = false,
            ) do t
                d_zt_part(k, z, t) * V(l, t) / t
            end

        res_az_thin_z_thin = res_az_thin_z_thin_part_1 + res_az_thin_z_thin_part_2

        if iswide(z)
            # Add enclosures of integral from az to az_thin and from
            # z_thin to z.
            # IMPROVE: Get better enclosures of this using mean value
            # theorem?
            (res_az_thin_z_thin) +
            ((az_thin - az) * d(k, z, az) * V(l, az) / az) +
            ((z - z_thin) * d(k, z, z) * V(l, z) / z)
        else
            res_az_thin_z_thin # We integrated everything
        end
    end

    return res_0_az + res_az_z
end

function integral_d_k_V_l_other_limit(k::Int, l::Int, z::Acb, b::Acb)
    a = 1e-5b

    # Integrate from 0 to a
    # TODO: Implement proper version of this
    if k == 3
        # TODO: Improve this enclosure using mean value theorem. That
        # requires evaluation with AcbSeries though.
        res_0_a_part1 = d_3_analytic_z(z, z) * V(l, z) / z
        res_0_a_part2 = d_3_analytic_conj_z(z, z) * V(l, z) / z

        res_0_a = Arblib.union(zero(a), a) * (res_0_a_part1 + res_0_a_part2)
    else
        res_0_a = Arblib.union(zero(a), a) * d(k, z, a) * V(l, a) / a
    end

    # Integrate from a to b
    # TODO: We need to verify analyticity for this to be correct,
    # which we might not have.
    res_a_b = if k == 3
        res_a_b_part1 = Arblib.integrate(
            a,
            b,
            atol = 1e-8,
            opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
            warn_on_no_convergence = false,
        ) do t
            d_3_zt_part(z, t) * V(l, t) / t
        end

        res_a_b_part2 =
            d_3_z_part(z) * Arblib.integrate(
                a,
                b,
                atol = 1e-8,
                opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
                warn_on_no_convergence = false,
            ) do t
                V(l, t) / t
            end

        res_a_b_part1 + res_a_b_part2
    else
        Arblib.integrate(
            a,
            b,
            atol = 1e-8,
            opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 0, 0),
            warn_on_no_convergence = false,
        ) do t
            d(k, z, t) * V(l, t) / t
        end
    end

    return res_0_a + res_a_b
end

function k_inv_N(inv_N)
    R = Arb("0.95")
    E_I = 5 // 2 * inv_N + 9 // 2 * inv_N^2 + 5inv_N^3 + 60inv_N^4
    λ = λ_disc()

    a₀ = _c_N(inv_N) / (sqrt(λ) * besselj1(sqrt(λ)))

    return (
        2Arb(π) * a₀^2 * R^2 / 2 *
        (besselj0(R * sqrt(λ_approx(inv_N)))^2 + besselj1(R * sqrt(λ_approx(inv_N)))^2) -
        4Arb(π) * a₀ * E_I * R^2 / 2 *
        hypgeom0f1_regularized(Arb(2), -Arb(1 // 4) * R^2 * λ_approx(inv_N)) - R^2 * E_I^2
    )
end

function ϵ_prime(inv_N)
    ϵ = 232inv_N^6 + 61inv_N^7 + 182inv_N^8

    sqrt(Arb(π)) * ϵ / sqrt(k_inv_N(inv_N))
end

function λ_sup_m_λ_inf(inv_N::Union{Arb,ArbSeries})
    # inv(N + 1)
    inv_Np1 = inv_N / (1 + inv_N)

    λ_disc() * (
        λ_approx_div_λ(inv_N) / (1 + ϵ_prime(inv_N)) -
        λ_approx_div_λ(inv_Np1) / (1 - ϵ_prime(inv_Np1))
    )
end

function q_sup_m_q_inf(inv_N::Union{Arb,ArbSeries})
    # inv(N + 1) and inv(N + 2)
    inv_Np1 = inv_N / (1 + inv_N)
    inv_Np2 = inv_Np1 / (1 + inv_Np1)

    λ_approx(inv_N) * (1 - ϵ_prime(inv_Np1)) / ((1 + ϵ_prime(inv_N)) * λ_approx(inv_Np1)) -
    λ_approx(inv_Np1) * (1 + ϵ_prime(inv_Np2)) /
    ((1 - ϵ_prime(inv_Np1)) * λ_approx(inv_Np2))
end
