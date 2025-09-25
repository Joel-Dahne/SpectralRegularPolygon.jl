_polylog(s::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}) = Arblib.polylog!(zero(z), s, z)

function polylog(s::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef})
    if s == 1
        return -log(1 - z)
    elseif Arblib.contains_zero(z) && 0.25 < abs_ubound(Arb, z) < 1
        # The Flint implementation fails when abs_ubound(Arb, z) is
        # too large. For that reason we use the direct bound when it
        # is a bit larger.
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(1, zᵤ) * zᵤ)
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
                throw(ArgumentError("only supports degrees up to 2"))
            end
        end
    end

    return ArbExtras.compose_zero!(res, res, z)
end


function S_integrand(n::Int, z::Arblib.AcbOrRef, t::Arblib.AcbOrRef)
    if n == 2
        # Note that n = 2 is the only case when the function is
        # bounded at t = 0.
        if Arblib.contains_zero(t)
            fx_div_x(Acb(t)) do t
                -4log(1 - t * z)
            end
        else
            return -4log(1 - t * z) / t
        end
    elseif 3 <= n <= 5
        logt = log(t)

        # We write the function as a polynomial in log(1 - t * z)
        p = AcbPoly()
        if n == 3
            p[1] = -12logt
            p[2] = 12
        elseif n == 4
            p[1] = -24logt^2
            p[2] = 48logt
            p[3] = -32
        elseif n == 5
            p[1] = -40logt^3
            p[2] = 120logt^2
            p[3] = -160logt
            p[4] = 80
        end

        return p(log(1 - t * z)) / t
    else
        ArbExtras.derivative_function(n) do inv_N
            inv_N * t^inv_N * ((1 - t * z)^-2inv_N - 1) / t
        end(Arb(0))
    end
end

function S(n::Int, z::Arblib.AcbOrRef)
    a = n == 2 ? Arb(0) : Arb(1e-8)
    b = Arblib.contains(z, Acb(1)) ? Arb(1) - 1e-8 : Arb(1)

    # Integrate from 0 to a
    res_0_a = if iszero(a)
        zero(z) # Nothing to integrate
    else
        let t = Arb((0, a))
            # Enclosure of log(1 - t * z) / t
            log_1mtz_div_t = fx_div_x(Acb(t)) do t
                log(1 - t * z)
            end
            log1mtz = log(1 - t * z)

            if n == 3
                (-12integral_log(1, a) + 12integral_log(0, a) * log1mtz) * log_1mtz_div_t / factorial(3)
            elseif n == 4
                (
                    -24integral_log(2, a) + 48integral_log(1, a) * log1mtz -
                    32integral_log(0, a) * log1mtz^2
                ) * log_1mtz_div_t / factorial(4)
            elseif n == 5
                (
                    -40integral_log(3, a) + 120integral_log(2, a) * log1mtz -
                    160integral_log(1, a) * log1mtz^2 +
                    80integral_log(0, a) * log1mtz^3
                ) * log_1mtz_div_t / factorial(5)
            else
                throw(ArgumentError("integration around zero not implemented for n > 5"))
            end
        end
    end

    # Integrate from a to b
    res_a_b =
        Arblib.integrate(
            t -> S_integrand(n, z, t),
            a,
            b,
            warn_on_no_convergence = false,
            atol = max(Arblib.radius(abs(res_0_a)), 2.0^-precision(z)),
            opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
        ) / factorial(n)

    # Integrate from b to 1
    res_b_1 = if isone(b)
        zero(z) # Nothing to integrate
    else
        let t = Arblib.union(b, Arb(1))
            if n == 2
                -2 / t * integral_log_1mtz(z, 1, b)

                -4integral_log_1mtz(z, 1, b) / (factorial(2) * t)
            elseif n == 3
                (-12log(t) * integral_log_1mtz(z, 1, b) + 12integral_log_1mtz(z, 2, b)) / (factorial(3) * t)
            elseif n == 4
                (
                    -24log(t)^2 * integral_log_1mtz(z, 1, b) +
                    48log(t) * integral_log_1mtz(z, 2, b) -
                    32integral_log_1mtz(z, 3, b)
                ) / (factorial(4) * t)
            elseif n == 5
                (
                    -40log(t)^3 * integral_log_1mtz(z, 1, b) +
                    120log(t)^2 * integral_log_1mtz(z, 2, b) -
                    160log(t) * integral_log_1mtz(z, 3, b) +
                    80integral_log_1mtz(z, 4, b)
                ) / (factorial(5) * t)
            else
                throw(ArgumentError("error bound not implemented for n > 5"))
            end
        end
    end

    return res_0_a + res_a_b + res_b_1
end

"""
    polylog_r_div_z_bound(r::Int, a::Arb)

For `0 < a < 1`, return `C` such that for any multiple polylogarithm
with `r` parameters (with are required to be positive integers) and
complex `z` with `abs(z) <= a`, the absolute value is bounded by `C *
z`.

IMPROVE: We can get better bounds for specific values of the
parameters, in particular when all values are 1. Might be worth it to
implement those specific bounds.
"""
function polylog_r_div_z_bound(r::Int, a::Arb)
    0 < a < 1 || return indeterminate(a)
    return 1 / (1 - a)^r
end


function polylog_1_1(z::Arblib.AcbOrRef)
    return log(1 - z)^2 / 2
end

function polylog_1_2(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(2, zᵤ) * zᵤ)
    else
        return log(1 - z)^2 * log(z) / 2 + log(1 - z) * polylog(2, 1 - z) -
               polylog(3, 1 - z) + zeta(Arb(3))
    end
end

function polylog_1_3(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(2, zᵤ) * zᵤ)
    else
        # NOTE: This is not the same formula as in the paper. It comes
        # from a recurrence relationship.
        return -polylog(4, 1 - z) + polylog(4, z) + polylog(4, inv(1 - 1 / z)) -
               polylog(3, z) * log(1 - z) + log(1 - z)^4 / factorial(4) -
               log(z) * log(1 - z)^3 / factorial(3) +
               zeta(Arb(2)) * log(1 - z)^2 / factorial(2) +
               zeta(Arb(3)) * log(1 - z) +
               zeta(Arb(4))
    end
end

function polylog_2_1(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(2, zᵤ) * zᵤ)
    else
        return -log(1 - z) / 6 * (Arb(π)^2 + 6polylog(2, 1 - z)) + 2polylog(3, 1 - z) -
               2zeta(Arb(3))
    end
end

function polylog_2_2(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(2, zᵤ) * zᵤ)
    else
        # NOTE: This is not quite the same version as in the paper, it
        # is slightly improved for better enclosures.
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

function polylog_3_1(z::Arblib.AcbOrRef)
    return -polylog(2, z)^2 / 2 - log(1 - z) * polylog(3, z)
end

function polylog_1_1_1(z::Arblib.AcbOrRef)
    return log(1 - z)^3 / 6
end

function polylog_1_1_2(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(3, zᵤ) * zᵤ)
    else
        return Arb(π)^4 / 90 - (1 // 6) * log(1 - z)^3 * log(z) -
               1 // 2 * log(1 - z)^2 * polylog(2, 1 - z) + log(1 - z) * polylog(3, 1 - z) -
               polylog(4, 1 - z)
    end
end

function polylog_1_2_1(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(3, zᵤ) * zᵤ)
    else
        return -Arb(π)^4 / 30 +
               1 // 2 * log(1 - z)^2 * polylog(2, 1 - z) +
               3 * polylog(4, 1 - z) - log(1 - z) * (2 * polylog(3, 1 - z) + zeta(Arb(3)))
    end
end

function polylog_2_1_1(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(3, zᵤ) * zᵤ)
    else
        return Arb(π)^4 / 30 +
               Arb(π)^2 / 12 * log(1 - z)^2 +
               log(1 - z) * polylog(3, 1 - z) - 3 * polylog(4, 1 - z) +
               2 * (-Acb(0, π) + log(z - 1)) * zeta(Arb(3))
    end
end

function polylog_1_1_1_1(z::Arblib.AcbOrRef)
    return log(1 - z)^4 / 24
end

function polylog_1_1_1_2(z::Arblib.AcbOrRef)
    if Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(4, zᵤ) * zᵤ)
    else
        return 1 // 24 * (
            log(1 - z)^4 * log(z) + 4log(1 - z)^3 * polylog(2, 1 - z) -
            12log(1 - z)^2 * polylog(3, 1 - z) + 24log(1 - z) * polylog(4, 1 - z) -
            24polylog(5, 1 - z)
        )
    end
end
