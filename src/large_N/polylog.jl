_polylog(s::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}) = Arblib.polylog!(zero(z), s, z)

"""
    _polylog_unitdisc(s::Int, z::Arblib.AcbOrRef)

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
function _polylog_unitdisc(s::Int, z::Arblib.AcbOrRef)
    N = 10000

    res = sum(1:(N-1)) do k
        z^k / Arb(k)^s
    end

    tail = abs(Arblib.polygamma!(zero(z), one(z), Acb(N)))

    return Arblib.add_error!(res, tail)
end

function polylog(s::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef})
    if Arblib.contains_zero(z) && 0.25 < abs_ubound(Arb, z) < 1
        # The Flint implementation fails when abs_ubound(Arb, z) is
        # too large. For that reason we use the direct bound when it
        # is a bit larger.
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^2 * zᵤ)
    elseif iswide(z) && !Arblib.contains_zero(z)
        # Explicit use of mean value theorem
        z_mid = Arblib.midpoint(Arblib._nonreftype(typeof(z)), z)

        f_mid = _polylog(s, z_mid)
        df = _polylog(s - 1, z) / z

        return Arblib.add_error(f_mid, abs(df) * abs(z - z_mid))
    else
        return _polylog(s, z)
    end
end

function polylog(s::Int, z::Union{ArbSeries,AcbSeries})
    z0 = z[0]
    res = zero(z)

    res[0] = polylog(s, z0)

    if length(z) > 1
        res[1] = polylog(s - 1, z0) / z0

        if length(z) > 2
            res[2] = (polylog(s - 2, z0) / z0 - res[1]) / 2z0

            if length(z) > 3
                error("not implemented")
            end
        end
    end

    return ArbExtras.compose_zero!(res, res, z)
end

"""
    polylog_unitdisc(s::Int, z::Arblib.AcbOrRef)

Compute `polylog(s, z)` assuming that `abs(z) <= 1`. If `z` overlaps
`1` it uses the (slowly converging) series expansion at zero through
`_polylog_unitdisc`, otherwise it falls back to `polylog`.
"""
function polylog_unitdisc(s::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef})
    if Arblib.contains(z, one(z))
        # TODO: We here assume that this is only ever called with
        # z in the unit disc.
        return _polylog_unitdisc(s, z)
    else
        return polylog(s, z)
    end
end

polylog_unitdisc(s::Int, z) = polylog(s, z)

lerch_phi(z::Acb, s::Int, a::Int) = Arblib.dirichlet_lerch_phi!(zero(z), z, Acb(s), Acb(a))

function S(n::Int, p::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef,ArbSeries,AcbSeries})
    s = n + 1
    if p == 1
        return polylog(s, z)
    elseif n == 1 && p == 2
        return polylog_1_2(z)
    elseif n == 1 && p == 3
        return polylog_1_1_2(z)
    elseif n == 1 && p == 4
        return polylog_1_1_1_2(z)
    elseif n == 2 && p == 2
        return polylog_1_3(z)
    elseif false #n == 2 && p == 3
        return polylog_1_1_3(z) # TODO: Implement this
    elseif false #n == 2 && p == 3
        return polylog_1_4(z) # TODO: Implement this
    elseif n == 1
        # TODO: We currently use hard coded versions above. It might
        # or might not be beneficial to use this recurrence.

        # Use recurrence relation from Proposition 2 in
        # https://arxiv.org/pdf/1908.04770. Using that S(1, p, 1) =
        # zeta(1 + p).
        # TODO: Handle z overlapping one
        # PROVE: The Proposition requires z != 1, but that doesn't seem required.
        (-1)^p // factorial(p) * log(z) * log(1 - z)^p + zeta(Arb(1 + p)) -
        sum(0:(p-1)) do k
            (-1)^k // factorial(k) * log(1 - z)^k * S(p - k, 1, 1 - z)
        end
    elseif p == 2
        # FIXME: This is not rigorous and converges slowly for abs(z)
        # = 1. Should prefer to rewrite in terms of polylog using
        # recurrence.
        sum(1:1000) do k₁
            z^(k₁ + 1) / k₁ * lerch_phi(z, s, k₁ + 1)
        end
    elseif p == 3
        # FIXME: Same as above
        sum(1:50) do k₁
            sum((k₁+1):51) do k₂
                z^(k₂ + 1) / k₂ * lerch_phi(z, s, k₂ + 1)
            end / k₁
        end
    elseif p == 4
        # FIXME: Same as above
        sum(1:20) do k₁
            sum((k₁+1):21) do k₂
                sum((k₂+1):22) do k₃
                    z^(k₃ + 1) / k₃ * lerch_phi(z, s, k₃ + 1)
                end / k₂
            end / k₁
        end
    end
end

function S_unitdisc(n::Int, p::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef})
    s = n + 1
    if p == 1
        return polylog_unitdisc(s, z)
    elseif n == 1 && p == 2
        return polylog_1_2(z)
    elseif n == 1 && p == 3
        return polylog_1_1_2(z)
    elseif n == 1 && p == 4
        return polylog_1_1_1_2(z)
    elseif n == 2 && p == 2
        return polylog_1_3(z)
    elseif false #n == 2 && p == 3
        return polylog_1_1_3(z) # TODO: Implement this
    elseif false #n == 2 && p == 3
        return polylog_1_4(z) # TODO: Implement this
    elseif n == 1
        # TODO: We currently use hard coded versions above. It might
        # or might not be beneficial to use this recurrence.

        # Use recurrence relation from Proposition 2 in
        # https://arxiv.org/pdf/1908.04770. Using that S(1, p, 1) =
        # zeta(1 + p).
        # TODO: Handle z overlapping one
        # PROVE: The Proposition requires z != 1, but that doesn't seem required.
        (-1)^p // factorial(p) * log(z) * log(1 - z)^p + zeta(Arb(1 + p)) -
        sum(0:(p-1)) do k
            (-1)^k // factorial(k) * log(1 - z)^k * S(p - k, 1, 1 - z)
        end
    elseif p == 2
        # FIXME: This is not rigorous and converges slowly for abs(z)
        # = 1. Should prefer to rewrite in terms of polylog using
        # recurrence.
        sum(1:1000) do k₁
            z^(k₁ + 1) / k₁ * lerch_phi(z, s, k₁ + 1)
        end
    elseif p == 3
        # FIXME: Same as above
        sum(1:50) do k₁
            sum((k₁+1):51) do k₂
                z^(k₂ + 1) / k₂ * lerch_phi(z, s, k₂ + 1)
            end / k₁
        end
    elseif p == 4
        # FIXME: Same as above
        sum(1:20) do k₁
            sum((k₁+1):21) do k₂
                sum((k₂+1):22) do k₃
                    z^(k₃ + 1) / k₃ * lerch_phi(z, s, k₃ + 1)
                end / k₂
            end / k₁
        end
    end
end

function S(n::Int, z)
    sum(1:(n-1)) do j
        (-1)^(j - 1) * 2^(n - j) * S(j, n - j, z)
    end
end

function S_integral(n::Int, z::Arblib.AcbOrRef)
    b = Arblib.contains(z, Acb(1)) ? Arb(1) - 1e-8 : Arb(1)

    # Integrate from 0 to b
    res_0_b = if n == 2
        # In this case the integrand is bounded at t = 0, so we can
        # integrate from 0 to b directly.
        Arblib.integrate(
            0,
            b,
            warn_on_no_convergence = false,
            opts = Arblib.calc_integrate_opt_struct(0, 1_000, 0, 0, 0),
        ) do t
            if Arblib.contains_zero(t)
                fx_div_x(Acb(t)) do t
                    -2log(1 - t * z)
                end
            else
                -2log(1 - t * z) / t
            end
        end
    else
        a = Arb(1e-8)

        # Integrate from a to b
        res_a_b =
            Arblib.integrate(
                a,
                b,
                warn_on_no_convergence = false,
                opts = Arblib.calc_integrate_opt_struct(0, 1_000, 0, 0, 0),
            ) do t
                ArbExtras.derivative_function(n) do inv_N
                    inv_N * t^inv_N * ((1 - t * z)^-2inv_N - 1) / t
                end(Arb(0))
            end / factorial(n)

        # Integrate from 0 to a
        res_0_a = let t = Arb((0, a))
            # Enclosure of log(1 - t * z) / t for t in [0, a]
            log_removable = fx_div_x(Acb(t)) do t
                log(1 - t * z)
            end

            if n == 3
                -2(integral_log(1, a) - integral_log(0, a) * log(1 - t * z)) * log_removable
            elseif n == 4
                (
                    -3integral_log(2, a) + 6integral_log(1, a) * log(1 - t * z) -
                    4integral_log(0, a) * log(1 - t * z)^2
                ) * log_removable / 3
            elseif n == 5
                (
                    -integral_log(3, a) + 3integral_log(2, a) * log(1 - t * z) -
                    4integral_log(1, a) * log(1 - t * z)^2 +
                    2integral_log(0, a) * log(1 - t * z)^3
                ) * log_removable / 3
            else
                throw(ArgumentError("error bound not implemented for n > 5"))
            end
        end

        res_0_a + res_a_b
    end

    # Integrate from b to 1
    res_b_1 = if isone(b)
        zero(res_0_b) # Nothing to integrate
    else
        let t = Arblib.union(b, Arb(1))
            if n == 2
                -2 / t * integral_log_1mtz(z, 1, b)
            elseif n == 3
                -2 / t * (log(t) * integral_log_1mtz(z, 1, b) - integral_log_1mtz(z, 2, b))
            elseif n == 4
                1 / 3t * (
                    -3log(t)^2 * integral_log_1mtz(z, 1, b) +
                    6log(t) * integral_log_1mtz(z, 2, b) - 4integral_log_1mtz(z, 3, b)
                )
            elseif n == 5
                1 / 3t * (
                    -log(t)^3 * integral_log_1mtz(z, 1, b) +
                    3log(t)^2 * integral_log_1mtz(z, 2, b) -
                    4log(t) * integral_log_1mtz(z, 3, b) + 2integral_log_1mtz(z, 4, b)
                )
            else
                throw(ArgumentError("error bound not implemented for n > 5"))
            end
        end
    end

    return res_0_b + res_b_1
end

# polylog(2, z) that allows evaluation around z = 1
polylog2(z::Acb) =
    if Arblib.contains(z, one(z))
        S_integral(2, z) / 2
    else
        polylog(2, z)
    end

function S_unitdisc(n::Int, z::Acb)
    sum(1:(n-1)) do j
        (-1)^(j - 1) * 2^(n - j) * S_unitdisc(j, n - j, z)
    end
end

function polylog_1_2(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^2 * zᵤ)
    elseif Arblib.contains(z, Acb(1))
        return indeterminate(z) # TODO: Implement this
    else
        return log(1 - z)^2 * log(z) / 2 + log(1 - z) * polylog(2, 1 - z) -
               polylog(3, 1 - z) + zeta(Arb(3))
    end
end

# NOTE: This is not the same version as in the paper
function polylog_1_3(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^2 * zᵤ)
    elseif Arblib.contains(z, Acb(1))
        return indeterminate(z) # TODO: Implement this
    else
        return -polylog(4, 1 - z) + polylog(4, z) + polylog(4, inv(1 - 1 / z)) -
               polylog(3, z) * log(1 - z) + log(1 - z)^4 / factorial(4) -
               log(z) * log(1 - z)^3 / factorial(3) +
               zeta(Arb(2)) * log(1 - z)^2 / factorial(2) +
               zeta(Arb(3)) * log(1 - z) +
               zeta(Arb(4))
    end
end

function polylog_1_3_v2(z::Arblib.AcbOrRef)
    return (1 // 360) * (
               Arb(π)^4 +
               15 * (
                   -6 * log(1 - z)^2 * log(z)^2 + 8 * log(1 - z) * log(z)^3 + log(z)^4 -
                   12 * log(1 - z) * log(z)^2 * log(-z / (-1 + z)) -
                   4 * log(z)^3 * log(-z / (-1 + z)) + 6 * log(z)^2 * log(-z / (-1 + z))^2 -
                   4 * log(1 / (1 - z)) * log(-z / (-1 + z))^3 -
                   4 * log(z) * log(-z / (-1 + z))^3 +
                   log(-z / (-1 + z))^4 +
                   12 * (log(z) - log(-z / (-1 + z)))^2 * polylog(2, 1 - z) -
                   12 * polylog(2, 1 - z)^2 +
                   12 *
                   log(z) *
                   (log(z) - 2 * (log(1 - z) + log(-z / (-1 + z)))) *
                   polylog(2, z) - 12 * log(-z / (-1 + z))^2 * polylog(2, z / (-1 + z)) -
                   48 * log(z) * polylog(3, 1 - z) +
                   24 * log(-z / (-1 + z)) * polylog(3, 1 - z) +
                   24 * log(1 - z) * polylog(3, z) +
                   24 * log(-z / (-1 + z)) * polylog(3, z) +
                   24 * log(-z / (-1 + z)) * polylog(3, z / (-1 + z)) +
                   24 * polylog(4, 1 - z) - 24 * polylog(4, z) -
                   24 * polylog(4, z / (-1 + z)) + 24 * log(z) * zeta(Arb(3))
               )
           ) - (11 // 360) * Arb(π)^4 +
           (1 // 12) * (
               -3 * log(1 / z)^4 - 4 * log(1 - z)^3 * (log(1 / z) - 2 * log(z)) -
               2 *
               log(1 - z) *
               (2 * Arb(π)^2 * log(1 / z) + 6 * log(1 / z)^3 + Arb(π)^2 * log(z)) -
               2 *
               log(1 - z)^2 *
               (Arb(π)^2 + 6 * log(1 / z)^2 - 6 * log(1 / z) * log(z) - 3 * log(z)^2)
           ) +
           (1 // 2) * polylog(2, 1 - z)^2 - log(-1 + 1 / z)^2 * polylog(2, (-1 + z) / z) +
           (log(-1 + 1 / z)^2 + log(1 - z) * log(z)) * polylog(2, z) +
           2 * (log(-1 + 1 / z) + log(z)) * polylog(3, 1 - z) +
           2 * log(-1 + 1 / z) * polylog(3, (-1 + z) / z) +
           2 * log(1 / z) * polylog(3, z) - 2 * polylog(4, 1 - z) -
           2 * polylog(4, (-1 + z) / z) + 2 * polylog(4, z)
end

function polylog_2_1(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^2 * zᵤ)
    elseif Arblib.contains(z, Acb(1))
        return indeterminate(z) # TODO: Implement this
    else
        return -log(1 - z) / 6 * (Arb(π)^2 + 6polylog(2, 1 - z)) + 2polylog(3, 1 - z) -
               2zeta(Arb(3))
    end
end

# NOTE: This is not quite the same version as in the paper, it is
# slightly improved for better enclosures.
function polylog_2_2(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^2 * zᵤ)
    elseif Arblib.contains(z, Acb(1))
        return indeterminate(z) # TODO: Implement this
    else
        return -Arb(π)^4 / 36 + polylog(2, 1 - z)^2 + (Arb(π)^2 / 6) * polylog(2, z) -
               2 * log(z) * zeta(Arb(3)) - (
            -(11 * Arb(π)^4 / 360) +
            (1 // 12) * (
                -3 * log(z)^4 - 4 * log(1 - z)^3 * (-3log(z)) -
                2 * log(1 - z) * (-Arb(π)^2 * log(z) - 6 * log(z)^3) -
                2 * log(1 - z)^2 * (Arb(π)^2 + 12 * log(z)^2 - 3 * log(z)^2)
            ) +
            1 // 2 * polylog(2, 1 - z)^2 - log(-1 + 1 / z)^2 * polylog(2, 1 - 1 / z) +
            (log(-1 + 1/z)^2 + log(1 - z) * log(z)) * polylog(2, z) +
            2 * log(-1 + 1/z) * polylog(3, 1 - z) +
            2 * log(-1 + 1/z) * polylog(3, 1 - 1 / z) +
            -2 * log(z) * polylog(3, z) - 2 * polylog(4, 1 - z) -
            2 * polylog(4, 1 - 1 / z) + 2 * polylog(4, z)
        )
    end
end

# NOTE: Not finite for z = 1
function polylog_3_1(z::Arblib.AcbOrRef)
    return -polylog(2, z)^2 / 2 - log(1 - z) * polylog(3, z)
end

function polylog_1_1_2(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^3 * zᵤ)
    elseif Arblib.contains(z, Acb(1))
        return indeterminate(z) # TODO: Implement this
    else
        return Arb(π)^4 / 90 - (1 // 6) * log(1 - z)^3 * log(z) -
               1 // 2 * log(1 - z)^2 * polylog(2, 1 - z) + log(1 - z) * polylog(3, 1 - z) -
               polylog(4, 1 - z)
    end
end

function polylog_1_2_1(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^3 * zᵤ)
    else
        return -Arb(π)^4 / 30 +
               1 // 2 * log(1 - z)^2 * polylog(2, 1 - z) +
               3 * polylog(4, 1 - z) - log(1 - z) * (2 * polylog(3, 1 - z) + zeta(Arb(3)))
    end
end

# NOTE: Not finite for z = 1
function polylog_2_1_1(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^3 * zᵤ)
    else
        return Arb(π)^4 / 30 +
               Arb(π)^2 / 12 * log(1 - z)^2 +
               log(1 - z) * polylog(3, 1 - z) - 3 * polylog(4, 1 - z) +
               2 * (-Acb(0, π) + log(z - 1)) * zeta(Arb(3))
    end
end

# NOTE: Not finite for z = 1
function polylog_1_1_1_1(z::Arblib.AcbOrRef)
    return (1 // 24) * log(1 - z)^4
end

function polylog_1_1_1_2(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        C = 1 / (1 - zᵤ)
        return add_error(zero(z), C^4 * zᵤ)
    elseif Arblib.contains(z, Acb(1))
        return indeterminate(z) # TODO: Implement this
    else
        return 1 // 24 * (
            log(1 - z)^4 * log(z) + 4log(1 - z)^3 * polylog(2, 1 - z) -
            12log(1 - z)^2 * polylog(3, 1 - z) + 24log(1 - z) * polylog(4, 1 - z) -
            24polylog(5, 1 - z)
        )
    end
end
