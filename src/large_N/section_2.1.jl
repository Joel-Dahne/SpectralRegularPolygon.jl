polylog(s::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}) = Arblib.polylog!(zero(z), s, z)

"""
    polylog_disc(s::Int, z::Arblib.AcbOrRef)

Compute `polylog(s, z)` assuming that `abs(z) <= 1`. It is intended to
be used for `z` overlapping `1`, otherwise there are much more
efficient methods.

It uses the power series expansion at `z = 0` and bounds the tail by
```
sum(k -> 1 / k^s, N:Inf) = polygamma(1, N)
```

IMPROVE: The bound is very slowly converging. It might be sufficient
for what we need though.
"""
function polylog_disc(s::Int, z::Arblib.AcbOrRef)
    N = 10000

    res = sum(1:(N-1)) do k
        z^k / Arb(k)^s
    end

    tail = abs(Arblib.polygamma!(zero(z), one(z), Acb(N)))

    return Arblib.add_error!(res, tail)
end

λ_disc() = ArbExtras.refine_root(besselj0, Arb((sqrt(Arf(5.7)), sqrt(Arf(5.8)))))^2

function S(n::Int, p::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef})
    if p == 1
        s = n + 1
        if Arblib.contains(z, one(z))
            # TODO: We here assume that this is only ever called with
            # z in the unit disc.
            return polylog_disc(s, z)
        else
            return polylog(s, z)
        end
    else
        #@warn "S(n, p) not implemented for p = $p"
        return zero(z)
    end
end

function S(n::Int, z)
    sum(1:(n-1)) do j
        (-1)^(j - 1) * 2^(n - j) * S(j, n - j, z)
    end
end

c_2(z) = (S(2, z) + S(2, conj(z))) / 2

c_3(z) = (S(3, z) + S(3, conj(z))) / 2

c_4(z) = -(S(2, z) - S(2, conj(z)))^2 / 8 + (S(4, z) + S(4, conj(z))) / 2

c_5(z) = (-(S(2, z) - S(2, conj(z))) + 2(S(5, z) + S(5, conj(z) - 2λ_disc() * zeta(5))))


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
               polylog(3, 1 - z),
               zeta(oftype(λ, 3)),
           ) +
           (λ - 4) * (
               -log(1 - z) / 6 * (oftype(λ, π)^2 + 6polylog(2, 1 - z)) +
               2polylog(3, 1 - z) - 2zeta(oftype(λ, 3))
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
           2^4 * polylog_1_1_1(z) +
           2λ * zeta(oftype(λ, 3)) * polylog(1, z)
end

function polylog_1_3(z)
    return Arb(π)^4 / 360 +
           (1 // 24) * (
        -90 * polylog(2, 1 - z)^2 + 90 * polylog(2, z)^2 -
        360 * log(z) * polylog(3, 1 - z) +
        180 * log(z / (z - 1)) * polylog(3, 1 - z) +
        180 * (log(1 - z) + log(z / (z - 1))) * polylog(3, z) +
        180 * log(z / (z - 1)) * polylog(3, z / (z - 1)) +
        360 * (log(z) - log(z / (z - 1)))^2 * polylog(2, 1 - z) +
        360 * log(z) * (log(z) - 2 * log(1 - z) - 2 * log(z / (z - 1))) * polylog(2, z) -
        360 * log(z / (z - 1))^2 * polylog(2, z / (z - 1)) + 360 * polylog(4, 1 - z) -
        360 * polylog(4, z) - 360 * polylog(4, z / (z - 1)) - 360 * log(z) * zeta(Arb(3)) +
        180 * log(z / (z - 1)) * zeta(Arb(3)) -
        90 * (
            -6 * log(1 - z)^2 * log(z)^2 + 8 * log(1 - z) * log(z)^3 + log(z)^4 -
            12 * log(1 - z) * log(z)^2 * log(z / (z - 1)) -
            4 * log(z)^3 * log(z / (z - 1)) + 6 * log(z)^2 * log(z / (z - 1))^2 -
            4 * log(1 - z) * log(z / (z - 1))^3 - 4 * log(z) * log(z / (z - 1))^3 +
            log(z / (z - 1))^4
        )
    )
end

function polylog_2_2(z)
    return -Arb(π)^4 / 36 +
           polylog(2, 1 - z)^2 +
           (Arb(π)^2 / 6) * polylog(2, z) +
           2 * log(z) * polylog(3, 1 - z) - 2 * log(z) * zeta(Arb(3)) - (
        (11 * Arb(π)^4 / 360) +
        (1 // 12) * (
            -3 * log(1/z)^4 - 4 * log(1 - z) * (log(1/z) - 2 * log(z)) -
            2 *
            log(1 - z) *
            (2 * Arb(π)^2 * log(1/z) + 6 * log(1/z)^3 + Arb(π)^2 * log(z)) -
            2 *
            log(1 - z)^2 *
            (Arb(π)^2 + 6 * log(1/z)^2 - 6 * log(1/z) * log(z) - 3 * log(z)^2)
        ) +
        1 // 2 * polylog(2, 1 - z)^2 - log(-1 + 1/z)^2 * polylog((z - 1)/z) +
        (log(-1 + 1/z)^2 + log(1 - z) * log(z)) * polylog(2, z) +
        2 * (log(-1 + 1/z) + log(z)) * polylog(3, 1 - z) +
        2 * log(-1 + 1/z) * polylog(3, (z - 1)/z) +
        2 * log(1/z) * polylog(3, z) - 2 * polylog(4, 1 - z) - 2 * polylog(4, (z - 1) / z) +
        2 * polylog(4, z)
    )
end

function polylog_3_1(z)
    return -polylog(2, z)^2 / 2 - log(1 - z) * polylog(3, z)
end

function polylog_2_1_1(z)
    return Arb(π)^4 / 30 + Arb(π)^2 / 12 * log(1 - z)^2 + log(1 - z) * polylog(3, 1 - z) -
           3 * polylog(4, 1 - z) + 2 * (-Acb(0, π) + log(z - 1)) * zeta(Arb(3))
end

function polylog_1_2_1(z)
    return -Arb(π)^4 / 30 +
           1 // 2 * log(1 - z)^2 * polylog(2, 1 - z) +
           3 * polylog(4, 1 - z) - log(1 - z) * (2 * polylog(3, 1 - z) + zeta(Arb(3)))
end

function polylog_1_1_2(z)
    return Arb(π)^4 / 90 - (1 // 6) * log(1 - z)^3 * log(z) -
           1 // 2 * log(1 - z)^2 * polylog(2, 1 - z) + log(1 - z) * polylog(3, 1 - z) -
           polylog(4, 1 - z)
end

function polylog_1_1_1_1(z)
    return (1 // 24) * log(1 - z)^4
end

function polylog_1_3_v2(z)
    return Arb(π)^4 / 360 +
           (1 // 24) * (
        -90 * polylog(2, 1 - z)^2 + 90 * polylog(2, z)^2 -
        360 * log(z) * polylog(3, 1 - z) +
        180 * log(z / (z - 1)) * polylog(3, 1 - z) +
        180 * (log(1 - z) + log(z / (z - 1))) * polylog(3, z) +
        180 * log(z / (z - 1)) * polylog(3, z / (z - 1)) +
        360 * (log(z) - log(z / (z - 1)))^2 * polylog(2, 1 - z) +
        360 * log(z) * (log(z) - 2 * log(1 - z) - 2 * log(z / (z - 1))) * polylog(2, z) -
        360 * log(z / (z - 1))^2 * polylog(2, z / (z - 1)) + 360 * polylog(4, 1 - z) -
        360 * polylog(4, z) - 360 * polylog(4, z / (z - 1)) - 360 * log(z) * zeta(Arb(3)) +
        180 * log(z / (z - 1)) * zeta(Arb(3)) -
        90 * (
            -6 * log(1 - z)^2 * log(z)^2 + 8 * log(1 - z) * log(z)^3 + log(z)^4 -
            12 * log(1 - z) * log(z)^2 * log(z / (z - 1)) -
            4 * log(z)^3 * log(z / (z - 1)) + 6 * log(z)^2 * log(z / (z - 1))^2 -
            4 * log(1 - z) * log(z / (z - 1))^3 - 4 * log(z) * log(z / (z - 1))^3 +
            log(z / (z - 1))^4
        )
    )
end

function polylog_2_2_v2(z)
    return -Arb(π)^4 / 36 +
           polylog(2, 1 - z)^2 +
           Arb(π)^2 / 6 * polylog(2, z) +
           2 * log(z) * polylog(3, 1 - z) - 2 * log(z) * zeta(Arb(3)) - (
        (11 * Arb(π)^4) // 360 +
        (1 // 12) * (
            -3 * log(1 / z)^4 - 4 * log(1 - z)^3 * (log(1 / z) - 2 * log(z)) -
            2 *
            log(1 - z) *
            (2 * Arb(π)^2 * log(1 / z) + 6 * log(1 / z)^3 + Arb(π)^2 * log(z)) -
            2 *
            log(1 - z)^2 *
            (Arb(π)^2 + 6 * log(1 / z)^2 - 6 * log(1 / z) * log(z) - 3 * log(z)^2)
        ) +
        (1 // 2) * polylog(2, 1 - z)^2 - log(-1 + 1 / z)^2 * polylog(2, (z - 1) / z) +
        (log(-1 + 1 / z)^2 + log(1 - z) * log(z)) * polylog(2, z) +
        2 * (log(-1 + 1 / z) + log(z)) * polylog(3, 1 - z) +
        2 * log(-1 + 1 / z) * polylog(3, (z - 1) / z) +
        2 * log(1 / z) * polylog(3, z) - 2 * polylog(4, 1 - z) -
        2 * polylog(4, (z - 1) / z) + 2 * polylog(4, z)
    )
end

function polylog_3_1_v2(z)
    return -(1 // 2) * polylog(2, z)^2 - log(1 - z) * polylog(3, z)
end

function polylog_2_1_1_v2(z)
    return Arb(π)^4 / 30 +
           (1 // 12) * Arb(π)^2 * log(1 - z)^2 +
           log(1 - z) * polylog(3, 1 - z) - 3 * polylog(4, 1 - z) +
           2 * (-im * Arb(π) + log(z - 1)) * zeta(Arb(3))
end

function polylog_1_2_1_v2(z)
    return -Arb(π)^4 / 30 +
           (1 // 2) * log(1 - z)^2 * polylog(2, 1 - z) +
           3 * polylog(4, 1 - z) - log(1 - z) * (2 * polylog(3, 1 - z) + zeta(Arb(3)))
end

function polylog_1_1_1_1_v2(z)
    return (1 // 24) * log(1 - z)^4
end

function polylog_1_3_v3(z)
    return Arb(π)^4 / 360 +
           (1 // 24) * (
        -90 * polylog(2, 1 - z)^2 + 90 * polylog(2, z)^2 -
        360 * log(z) * polylog(3, 1 - z) +
        180 * log(z / (z - 1)) * polylog(3, 1 - z) +
        180 * (log(1 - z) + log(z / (z - 1))) * polylog(3, z) +
        180 * log(z / (z - 1)) * polylog(3, z / (z - 1)) +
        360 * (log(z) - log(z / (z - 1)))^2 * polylog(2, 1 - z) +
        360 * log(z) * (log(z) - 2 * log(1 - z) - 2 * log(z / (z - 1))) * polylog(2, z) -
        360 * log(z / (z - 1))^2 * polylog(2, z / (z - 1)) + 360 * polylog(4, 1 - z) -
        360 * polylog(4, z) - 360 * polylog(4, z / (z - 1)) - 360 * log(z) * zeta(Arb(3)) +
        180 * log(z / (z - 1)) * zeta(Arb(3)) -
        90 * (
            -6 * log(1 - z)^2 * log(z)^2 + 8 * log(1 - z) * log(z)^3 + log(z)^4 -
            12 * log(1 - z) * log(z)^2 * log(z / (z - 1)) -
            4 * log(z)^3 * log(z / (z - 1)) + 6 * log(z)^2 * log(z / (z - 1))^2 -
            4 * log(1 - z) * log(z / (z - 1))^3 - 4 * log(z) * log(z / (z - 1))^3 +
            log(z / (z - 1))^4
        )
    )
end

function polylog_2_2_v3(z)
    return -Arb(π)^4 / 36 +
           polylog(2, 1 - z)^2 +
           (Arb(π)^2 / 6) * polylog(2, z) +
           2 * log(z) * polylog(3, 1 - z) - 2 * log(z) * zeta(Arb(3)) - (
               (11 * Arb(π)^4) / 360 +
               (1 // 12) * (
                   -3 * log(1 / z)^4 - 4 * log(1 - z)^3 * (log(1 / z) - 2 * log(z)) -
                   2 *
                   log(1 - z) *
                   (2 * Arb(π)^2 * log(1 / z) + 6 * log(1 / z)^3 + Arb(π)^2 * log(z)) -
                   2 *
                   log(1 - z)^2 *
                   (Arb(π)^2 + 6 * log(1 / z)^2 - 6 * log(1 / z) * log(z) - 3 * log(z)^2)
               )
           ) + (1 // 2) * polylog(2, 1 - z)^2 -
           log(-1 + 1 / z)^2 * polylog(2, (z - 1) / z) +
           (log(-1 + 1 / z)^2 + log(1 - z) * log(z)) * polylog(2, z) +
           2 * (log(-1 + 1 / z) + log(z)) * polylog(3, 1 - z) +
           2 * log(-1 + 1 / z) * polylog(3, (z - 1) / z) +
           2 * log(1 / z) * polylog(3, z) - 2 * polylog(4, 1 - z) -
           2 * polylog(4, (z - 1) / z) + 2 * polylog(4, z)
end

function polylog_3_1_v3(z)
    return -1 // 2 * polylog(2, z)^2 - log(1 - z) * polylog(3, z)
end

function polylog_2_1_1_v3(z)
    return Arb(π)^4 / 30 +
           (1 // 12) * Arb(π)^2 * log(1 - z)^2 +
           log(1 - z) * polylog(3, 1 - z) - 3 * polylog(4, 1 - z) +
           2 * (-im * Arb(π) + log(z - 1)) * zeta(Arb(3))
end

function polylog_1_2_1_v3(z)
    return -Arb(π)^4 / 30 +
           (1 // 2) * log(1 - z)^2 * polylog(2, 1 - z) +
           3 * polylog(4, 1 - z) - log(1 - z) * (2 * polylog(3, 1 - z) + zeta(Arb(3)))
end

function polylog_1_1_2_v3(z)
    return Arb(π)^4 / 90 - (1 // 6) * log(1 - z)^3 * log(z) -
           (1 // 2) * log(1 - z)^2 * polylog(2, 1 - z) + log(1 - z) * polylog(3, 1 - z) -
           polylog(4, 1 - z)
end

function polylog_1_1_1_1_v3(z)
    return (1 // 24) * log(1 - z)^4
end

c(N::Int) = sqrt(
    gamma(1 - Arb(1 // N))^2 * gamma(1 + Arb(2 // N)) /
    (gamma(1 + Arb(1 // N))^2 * gamma(1 - Arb(2 // N))),
)

c_inv_N(inv_N) =
    sqrt(gamma(1 - inv_N)^2 * gamma(1 + 2inv_N) / (gamma(1 + inv_N)^2 * gamma(1 - 2inv_N)))
