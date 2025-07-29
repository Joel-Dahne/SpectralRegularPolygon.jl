λ_disc() = ArbExtras.refine_root(besselj0, Arb((sqrt(Arf(5.7)), sqrt(Arf(5.8)))))^2
function λ_approx(inv_N::Arb)
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
    a = Acb(1e-8) # FIXME

    return Arblib.integrate(a, 1) do t
        t^inv_N * ((1 - t * z)^-2inv_N - 1) / t
    end
end

function F_N(inv_N::Arb, z::Acb)
    return 1 + inv_N * F_N_sub_1_mul_N(inv_N, z)
end

# c_2(z) = S(2, z) + S(2, conj(z))) / 2
# NOTE: This assumes that abs(z) == 1
function c_2(z::Acb)
    if Arblib.contains(z, Acb(1, 0))
        # Input overlaps branch cut
        return real(S_unitdisc(2, z))
    else
        return real(mean_value_theorem_bound(z -> S(2, z), z))
    end
end

# c_3(z) = (S(3, z) + S(3, conj(z))) / 2
function c_3(z::Acb)
    return real(mean_value_theorem_bound(z -> S(3, z), z))
end

# c_4(z) = -(S(2, z) - S(2, conj(z)))^2 / 8 + (S(4, z) + S(4, conj(z))) / 2
function c_4(z::Acb)
    S2_imag = imag(mean_value_theorem_bound(z -> S(2, z), z))
    S4_real = real(mean_value_theorem_bound(z -> S(4, z), z))

    return S2_imag^2 / 2 + S4_real
end

#c_5(z) =
#    (
#        -(S(2, z) - S(2, conj(z))) * (S(3, z) - S(3, conj(z))) +
#        2(S(5, z) + S(5, conj(z)) - 2λ_disc() * zeta(Arb(5)))
#    ) / 4
function c_5(z::Acb)
    S2_imag = imag(mean_value_theorem_bound(z -> S(2, z), z))
    S3_imag = imag(mean_value_theorem_bound(z -> S(3, z), z))
    S5_real = real(mean_value_theorem_bound(z -> S(5, z), z))

    return S2_imag * S3_imag + S5_real - λ_disc() * zeta(Arb(5))
end
_c_5_real(z) = S(5, z) - λ_disc() * zeta(Arb(5))
_c_5_imag(z) = imag(S(2, z)) * imag(S(3, z))

function T_2_bound(N₀::Int)
    # Compute all parts not depending on z
    inv_N = Arb((0, 1 // N₀))

    # (c_N - 1) * N^2
    factor_part1 = c_N_remainder_mul_N3(N₀) * inv_N

    # (4zeta(3) / N + (12 - 2λ) * zeta(5) / N^3) * c_N
    # IMPROVE: We can move some of the N to c_N to get better
    # enclosures.
    factor_part2 =
        (4zeta(Arb(3)) * inv_N + (12 - 2λ_disc()) * zeta(Arb(5)) * inv_N^3) * c_N(N₀)

    factor = factor_part1 + factor_part2

    return z -> let
        # Bound of main term
        part1 = let
            a = Acb(1e-8) # FIXME

            Arblib.integrate(
                a,
                1,
                atol = 1e-4,
                opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
            ) do t
                fx_div_x(Acb(inv_N), extra_degree = 2) do inv_N
                    t^inv_N * ((1 - t * z)^-2inv_N - 1)
                end / t
                #let inv_N = Arb(1 // 26)
                #    t^inv_N * ((1 - t * z)^-2inv_N - 1) / t / inv_N
                #end
            end
        end

        # Enclosure of other term
        part2 = factor * F_N(inv_N, z)

        return abs(part1) + abs(part2)
    end
end

# Simplified version of T_2. It is just (F_N(z) - 1) * N^2
function T_2_simplified(inv_N::Arb, z::Acb)
    #return F_N_sub_1_mul_N(inv_N, z) / inv_N

    # FIXME
    a = Acb(1e-8)

    #res1 = 2Arblib.integrate(a, b) do t
    #    -log(1 - t * z) / t
    #end
    res1 = 2polylog(2, z)

    res2 = if Arblib.contains_zero(inv_N)
        Arblib.integrate(
            a,
            1,
            atol = 1e-2Arblib.ubound(abs(res1)),
            opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
        ) do t
            fx_div_x(Acb(inv_N), 2, extra_degree = 2, force = true) do inv_N
                (t^inv_N * ((1 - t * z)^-2inv_N - 1) + 2log(1 - t * z) * inv_N)
            end / t
        end * inv_N
    else
        Arblib.integrate(a, 1) do t
            (t^inv_N * ((1 - t * z)^-2inv_N - 1) + 2log(1 - t * z) * inv_N) / t
        end / inv_N
    end

    return res1 + res2
end

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
        # (abs(F_N(z) - 1 - c_2(z) / N^2 - c_3(z) / N^3 - c_4(z) / N^4 - XXX)) * N^6
        term1 = zero(Arb)#indeterminate(Arb)

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

function g_dw3(w)
    λ = λ_disc()
    return λ^(3 // 2) * (3besselj1(w * sqrt(λ) - besselj(oftype(λ, 3), w * sqrt(λ)))) / 4
end


V_1(z) = -2log(1 - z)

V_2(z) = (λ_disc() / 2 - 2) * polylog(2, z) + 2log(1 - z)^2

function V_3(z)
    λ = λ_disc()
    return (λ^2 / 16 - λ + 2) * polylog(3, z) +
           (3λ - 12) * (
               log(1 - z)^2 * log(z) / 2 + log(1 - z) * polylog(2, 1 - z) -
               polylog(3, 1 - z) + zeta(Arb(3))
           ) +
           (λ - 4) * (
               -log(1 - z) / 6 * (Arb(π)^2 + 6polylog(2, 1 - z)) + 2polylog(3, 1 - z) -
               2zeta(Arb(3))
           ) - 4log(1 - z)^3 / 3
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
           2^4 * polylog_1_1_1_1(z) +
           2λ * zeta(Arb(3)) * polylog(1, z)
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

function d_0(z, t)
    return 1
end

function d_1(z, t)
    λ = λ_disc()
    return λ * log(t / z) / 4
end

function d_2(z, t)
    λ = λ_disc()
    return (1 // 64) * λ^2 * log(t / z)^2 +
           (1 // 8) * λ * (log(t / z)^2 + 4 * polylog(2, t) - 4polylog_unsafe(2, z))
end

function d_3(z, t)
    λ = λ_disc()

    # NOTE: This is an alternative formulation that is SLIGHTLY better
    # than the version in the paper.
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + (1 // 16) * λ + 1 // 6) * log(t / z)^2 +
            ((1 // 8) * λ + 1) * S(2, t) - (1 // 8) * λ * S(2, z) + S(2, conj(z))
        ) + S(3, t) - S(3, z)
    )

    # This is the version in the paper
    #return λ^3 * log(t / z)^3 / 2304 +
    #       (1 // 64) * λ^2 * log(t / z) * (log(t / z)^2 + 2S(2, t) - 2S(2, z)) +
    #       (1 // 24) *
    #       λ *
    #       (
    #           log(t / z)^3 + 6log(t / z) * S(2, t) + 6log(t / z) * S(2, conj(z)) + 6S(3, t) -
    #           6S(3, z)
    #       )
end

# Part of d_3 depending on both z and t
function d_3_p1(z, t)
    λ = λ_disc()
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + (1 // 16) * λ + 1 // 6) * log(t / z)^2 +
            ((1 // 8) * λ + 1) * S(2, t) - (1 // 8) * λ * S(2, z) + S(2, conj(z))
        ) + S(3, t)
    )
end

# Part of d_3 depending only on z
function d_3_p2(z)
    λ = λ_disc()
    return -λ / 4 * S(3, z)
end

# Part of d_3 analytic in z
function d_3_part_1(z, t)
    λ = λ_disc()
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + (1 // 16) * λ + 1 // 6) * log(t / z)^2 +
            ((1 // 8) * λ + 1) * S(2, t) - (1 // 8) * λ * S(2, z)
        ) + S(3, t) - S(3, z)
    )
end

# Part of d_3 analytic in conj(z)
function d_3_part_2(z, t)
    λ = λ_disc()
    return λ / 4 * log(t / z) * S(2, conj(z))
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

function integral_d_k_V_l(k::Int, l::Int, z::Acb)
    a = 1e-5Arblib.midpoint(Acb, z)
    b = Arblib.midpoint(Acb, z)

    # Integrate from 0 to a
    # TODO: Implement proper version of this
    if k == 3
        res1_part1 = mean_value_theorem_bound(z) do z
            d_3_part_1(z, z) * V(l, z) / z
        end
        res1_part2 = d_3_part_2(z, z) * V(l, z) / z

        res1 = Arblib.union(zero(a), a) * (res1_part1 + res1_part2)
    else
        res1 = Arblib.union(zero(a), a) * d(k, z, a) * V(l, a) / a
    end

    # Integrate from a to b
    # TODO: We need to verify analyticity for this to be correct,
    # which we might not have.
    if k == 3
        res2_part1 = Arblib.integrate(
            a,
            b,
            atol = 1e-8,
            opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
        ) do t
            d_3_p1(z, t) * V(l, t) / t
        end

        res2_part2 =
            mean_value_theorem_bound(d_3_p2, z) * Arblib.integrate(
                a,
                b,
                atol = 1e-8,
                opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
            ) do t
                V(l, t) / t
            end

        res2 = res2_part1 + res2_part2
    else
        res2 = Arblib.integrate(
            a,
            b,
            atol = 1e-8,
            opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
        ) do t
            d(k, z, t) * V(l, t) / t
        end
    end

    # Integrate from b to z by enclosing integrand on z and
    # multiplying by radius.
    # TODO: Verify that this is correct
    if k == 3
        res3_part1 = abs(z - b) * mean_value_theorem_bound(z) do z
            d_3_part_1(z, z) * V(l, z) / z
        end
        res3_part2 = d_3_part_2(z, z) * V(l, z) / z

        res3 = res3_part1 + res3_part2
    else
        res3 = abs(z - b) * mean_value_theorem_bound(z) do z
            d(k, z, z) * V(l, z) / z
        end
    end

    return res1 + res2 + res3
end

function integral_d_k_V_l_other_limit(k::Int, l::Int, z::Acb, b::Acb)
    a = 1e-5b

    # Integrate from 0 to a
    # TODO: Implement proper version of this
    if k == 3
        res1_part1 = mean_value_theorem_bound(z) do z
            d_3_part_1(z, z) * V(l, z) / z
        end
        res1_part2 = d_3_part_2(z, z) * V(l, z) / z

        res1 = Arblib.union(zero(a), a) * (res1_part1 + res1_part2)
    else
        res1 = Arblib.union(zero(a), a) * d(k, z, a) * V(l, a) / a
    end

    # Integrate from a to b
    # TODO: We need to verify analyticity for this to be correct,
    # which we might not have.
    if k == 3
        res2_part1 = Arblib.integrate(
            a,
            b,
            atol = 1e-8,
            opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
        ) do t
            d_3_p1(z, t) * V(l, t) / t
        end

        res2_part2 =
            mean_value_theorem_bound(d_3_p2, z) * Arblib.integrate(
                a,
                b,
                atol = 1e-8,
                opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
            ) do t
                V(l, t) / t
            end

        res2 = res2_part1 + res2_part2
    else
        res2 = Arblib.integrate(
            a,
            b,
            atol = 1e-8,
            opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 0, 0),
        ) do t
            d(k, z, t) * V(l, t) / t
        end
    end

    return res1 + res2
end

c(N::Int) = sqrt(
    gamma(1 - Arb(1 // N))^2 * gamma(1 + Arb(2 // N)) /
    (gamma(1 + Arb(1 // N))^2 * gamma(1 - Arb(2 // N))),
)

c_inv_N(inv_N) =
    sqrt(gamma(1 - inv_N)^2 * gamma(1 + 2inv_N) / (gamma(1 + inv_N)^2 * gamma(1 - 2inv_N)))
