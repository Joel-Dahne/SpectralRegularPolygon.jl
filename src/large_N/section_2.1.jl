λ_disc() = ArbExtras.refine_root(besselj0, Arb((sqrt(Arf(5.7)), sqrt(Arf(5.8)))))^2

# c_2(z) = S(2, z) + S(2, conj(z))) / 2
function c_2(z::Acb)
    return real(mean_value_theorem_bound(z -> S(2, z), z))
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
    #z = Arblib.midpoint(Acb, z)

    return λ^3 * log(t / z)^3 / 2304 +
           (1 // 64) * λ^2 * log(t / z) * (log(t / z)^2 + 2S(2, t) - 2S(2, z)) +
           (1 // 24) *
           λ *
           (
               log(t / z)^3 +
               6log(t / z) * S(2, t) +
               6log(t / z) * S(2, conj(z)) +
               6S(3, t) - 6S(3, z)
           )
end

function d_3_p1(z, t)
    λ = λ_disc()
    return λ^3 * log(t / z)^3 / 2304 +
           (1 // 64) * λ^2 * log(t / z) * (log(t / z)^2 + 2S(2, t) - 2S(2, z)) +
           (1 // 24) *
           λ *
           (log(t / z)^3 + 6log(t / z) * S(2, t) + 6log(t / z) * S(2, conj(z)) + 6S(3, t))
end

function d_3_p2(z)
    λ = λ_disc()
    return (1 // 24) * λ * (-6S(3, z))
end

function d_3_part_1(z, t)
    λ = λ_disc()
    return λ^3 * log(t / z)^3 / 2304 +
           (1 // 64) * λ^2 * log(t / z) * (log(t / z)^2 + 2S(2, t) - 2S(2, z)) +
           (1 // 24) * λ * (log(t / z)^3 + 6log(t / z) * S(2, t) + 6S(3, t) - 6S(3, z))
end

function d_3_part_2(z, t)
    λ = λ_disc()
    return (1 // 24) * λ * 6log(t / z) * S(2, conj(z))
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
    a = 0.001Arblib.midpoint(Acb, z)
    b = Arblib.midpoint(Acb, z)

    # Integrate from 0 to a
    # TODO: Implement proper version of this
    res1 = Arblib.union(zero(a), a) * d(k, z, a) * V(l, a) / a

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
            opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 1, 0),
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
    a = 0.001b

    # Integrate from 0 to a
    # TODO: Implement this
    res1 = zero(z)

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
            opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
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
