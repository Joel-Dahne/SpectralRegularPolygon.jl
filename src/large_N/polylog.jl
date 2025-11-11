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
        # For wide values of z the Flint implementation gives very
        # poor bounds. We make use of the mean value theorem to get
        # improved enclosures in this case.

        z_mid = midpoint(Arblib._nonreftype(typeof(z)), z)

        # Evaluate at the midpoint
        f_mid = _polylog(s, z_mid)

        # Bound derivative on the full interval
        df = _polylog(s - 1, z) / z

        return add_error(f_mid, abs(df) * abs(z - z_mid))
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

"""
    S_integrand(n::Int, z::Arblib.AcbOrRef, t::Arblib.AcbOrRef; analytic::Bool)

Compute the `n`th derivative of the function
```
inv_N * t^inv_N * ((1 - t * z)^-2inv_N - 1) / t
```
in terms of `inv_N`, evaluated at `inv_N = 0`. This is the integrand
in the integral formula for [`S`](@ref).

The derivative could be computed automatically using Taylor
expansions. For performance reasons we however evaluate it using
explicit expressions implemented for `2 <= n <= 5`.

If `analytic` is set, then the function returns an indeterminate
result if the function is not analytic on the entire interval `t`.
This follows the convention of [`Arblib.integrate`](@ref). For this we
have to check if `t` overlaps with a branch cut of the function (poles
automatically give indeterminate values). The branch cuts are given by
`t` or `1 - t * z` overlapping the interval ``(-∞, 0]``. In the
special case that `n = 2` there is not branch cut for `t` on this
interval.
"""
function S_integrand(n::Int, z::Arblib.AcbOrRef, t::Arblib.AcbOrRef; analytic::Bool)
    if analytic
        # Check if t overlaps the branch cut. The branch cut is for
        # either t or 1 - tz lying on the negative real axis. For the
        # special case that n = 2 there is no branch cut for t on the
        # negative real axis.

        # If n != 2, check if t overlaps the non-positive real axis.
        if n != 2 &&
           Arblib.contains_nonpositive(Arblib.realref(t)) &&
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
        throw(ArgumentError("only supports 2 <= n <= 5"))
    end
end

"""
    S(n::Int, z::Arblib.AcbOrRef)

Compute the function ``S_n(z)`` from the paper.

It is computed using the integral expression
```
S(n, z) = ∫ S_integral(n, z, t) dt
```
where the integral is taken from `0` to `1` and we refer to
[`S_integral(n, z, t)`](ref) for the precise integrand.

The integral is computed by splitting it into three parts, given by
integration on the intervals ``[0, a]``, ``[a, b]`` and ``[b, 1]``.
Note that in some cases `a` and `b` can be taken to be `0` and `1`
respectively. The interval ``[a, b]`` is handled by integrating
directly using [`Arblib.integrate`](@ref). The other two intervals are
handled by factoring out the bounded parts of the integrand and
explicitly integrating the rest.
"""
function S(n::Int, z::Arblib.AcbOrRef)
    2 <= n <= 5 || throw(ArgumentError("only supports 2 <= n <= 5"))

    a = n == 2 ? Arb(0) : Arb(1e-8)
    b = Arblib.contains(z, Acb(1)) ? Arb(1 - 1e-8) : Arb(1)

    # Note that res is S(n, z) * factorial(n), we divide away
    # factorial(n) at the end.

    # Integrate from 0 to a
    res_0_a = if iszero(a)
        zero(z) # Nothing to integrate
    else
        # Factor out bounds for log(1 - t * z) / t as well as log(1 -
        # t * z) from the explicit integrands. The resulting integrals
        # can then be computed explicitly using integral_log.

        # Enclosures of log(1 - t * z) / t and log(1 - t * z)
        log_1mtz_div_t, log1mtz = let t = Arb((0, a))
            fx_div_x(t -> log(1 - t * z), Acb(t)), log(1 - t * z)
        end

        if n == 3
            (-12integral_log(1, a) + 12integral_log(0, a) * log1mtz) * log_1mtz_div_t
        elseif n == 4
            (
                -24integral_log(2, a) + 48integral_log(1, a) * log1mtz -
                32integral_log(0, a) * log1mtz^2
            ) * log_1mtz_div_t
        elseif n == 5
            (
                -40integral_log(3, a) + 120integral_log(2, a) * log1mtz -
                160integral_log(1, a) * log1mtz^2 + 80integral_log(0, a) * log1mtz^3
            ) * log_1mtz_div_t
        end
    end

    # Integrate from a to b
    res_a_b = Arblib.integrate(
        (t; analytic) -> S_integrand(n, z, t; analytic),
        a,
        b,
        check_analytic = true,
        warn_on_no_convergence = false,
        atol = max(radius(abs(res_0_a)), 2.0^-precision(z)),
        opts = Arblib.calc_integrate_opt_struct(0, 2_000, 0, 0, 0),
    )

    # Integrate from b to 1
    res_b_1 = if isone(b)
        zero(z) # Nothing to integrate
    else
        # Factor out bounds for log(t) as well as 1 / t from the
        # explicit integrands. The resulting integrals can then be
        # computed explicitly using integral_log_1mtz.

        # Verify that the real and imaginary parts of log(1 - t * z)^m
        # don't change sign on the interval, for 1 <= m <= n - 1.
        let t = Arb((b, 1))
            C = Arblib.abs_ubound(Arb, 1 - t * z)
            C < 1 || return indeterminate(z)

            # We are now ensured that log(1 - t * z) lies in the strip
            # (-∞, log(C)) × im * [0, π] or (-∞, log(C)) × im * [π, 0]

            # The strip lies inside an angular sector with angle
            θ = atan(Arb(π), abs(log(C)))

            # To ensure that log(1 - t * z)^m lies inside a single
            # quadrant we have to verify that the angular sector which
            # it is contained in doesn't have an angle greater than π
            # / 2.
            (n - 1) * θ <= Arb(π) / 2 || return indeterminate(z)
        end

        # Enclosures of log(t) and 1 / t
        logt, invt = let t = Arb((b, 1))
            log(t), inv(t)
        end

        if n == 2
            -4integral_log_1mtz(z, 1, b) * invt
        elseif n == 3
            (-12logt * integral_log_1mtz(z, 1, b) + 12integral_log_1mtz(z, 2, b)) * invt
        elseif n == 4
            (
                -24logt^2 * integral_log_1mtz(z, 1, b) +
                48logt * integral_log_1mtz(z, 2, b) - 32integral_log_1mtz(z, 3, b)
            ) * invt
        elseif n == 5
            (
                -40logt^3 * integral_log_1mtz(z, 1, b) +
                120logt^2 * integral_log_1mtz(z, 2, b) -
                160logt * integral_log_1mtz(z, 3, b) + 80integral_log_1mtz(z, 4, b)
            ) * invt
        end
    end

    return (res_0_a + res_a_b + res_b_1) / factorial(n)
end

"""
    polylog_r_div_z_bound(r::Int, a::Arb)

For `0 < a < 1`, return `C` such that for any multiple polylogarithm
with `r` parameters (with are required to be positive integers) and
complex `z` with `abs(z) <= a`, the absolute value is bounded by `C *
z`.
"""
function polylog_r_div_z_bound(r::Int, a::Arb)
    0 < a < 1 || return indeterminate(a)
    return 1 / (1 - a)^r
end

# The expressions for the multiple polylogarithms below are given in
# the Appendix of the paper.

# All of these functions accept AcbSeries as input. This is only used
# in the tests to check if they expressions are correct by check that
# they satisfy the expected rules for derivatives.

function polylog_1_1(z::Union{Arblib.AcbOrRef,AcbSeries})
    return log(1 - z)^2 / 2
end

function polylog_1_2(z::Union{Arblib.AcbOrRef,AcbSeries})
    if z isa Arblib.AcbOrRef && Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(2, zᵤ) * zᵤ)
    else
        return 1 // 2 * log(1 - z)^2 * log(z) + log(1 - z) * polylog(2, 1 - z) -
               polylog(3, 1 - z) + zeta(Arb(3))
    end
end

function polylog_1_3(z::Union{Arblib.AcbOrRef,AcbSeries})
    if z isa Arblib.AcbOrRef && Arblib.contains_zero(z) && abs(z) < 1
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

function polylog_2_1(z::Union{Arblib.AcbOrRef,AcbSeries})
    if z isa Arblib.AcbOrRef && Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(2, zᵤ) * zᵤ)
    else
        return -1 // 6 * log(1 - z) * (Arb(π)^2 + 6polylog(2, 1 - z)) + 2polylog(3, 1 - z) -
               2zeta(Arb(3))
    end
end

function polylog_2_2(z::Union{Arblib.AcbOrRef,AcbSeries})
    if z isa Arblib.AcbOrRef && Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(2, zᵤ) * zᵤ)
    else
        return -log(1 - z)^3 * log(z) + 1 // 6 * log(1 - z)^2 * (Arb(π)^2 + 9 * log(z)^2) -
               1 // 6 * log(1 - z) * (Arb(π)^2 * log(z) + 6log(z)^3) +
               1 // 4 * log(z)^4 +
               -1 // 2 * polylog(2, 1 - z)^2 +
               log(-1 + 1 / z)^2 * polylog(2, 1 - 1 / z) -
               (log(-1 + 1/z)^2 + log(1 - z) * log(z)) * polylog(2, z) -
               2log(-1 + 1/z) * polylog(3, 1 - z) - 2log(-1 + 1/z) * polylog(3, 1 - 1 / z) -
               -2log(z) * polylog(3, z) +
               2polylog(4, 1 - z) +
               2polylog(4, 1 - 1 / z) - 2polylog(4, z) +
               polylog(2, 1 - z)^2 +
               Arb(π)^2 / 6 * polylog(2, z) - 2log(z) * zeta(Arb(3)) + Arb(π)^4 / 360
    end
end

function polylog_3_1(z::Union{Arblib.AcbOrRef,AcbSeries})
    return -log(1 - z) * polylog(3, z) - 1 // 2 * polylog(2, z)^2
end

function polylog_1_1_1(z::Union{Arblib.AcbOrRef,AcbSeries})
    return -log(1 - z)^3 / 6
end

function polylog_1_1_2(z::Union{Arblib.AcbOrRef,AcbSeries})
    if z isa Arblib.AcbOrRef && Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(3, zᵤ) * zᵤ)
    else
        return -1 // 6 * log(1 - z)^3 * log(z) - 1 // 2 * log(1 - z)^2 * polylog(2, 1 - z) +
               log(1 - z) * polylog(3, 1 - z) - polylog(4, 1 - z) + Arb(π)^4 / 90
    end
end

function polylog_1_2_1(z::Union{Arblib.AcbOrRef,AcbSeries})
    if z isa Arblib.AcbOrRef && Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(3, zᵤ) * zᵤ)
    else
        return 1 // 2 * log(1 - z)^2 * polylog(2, 1 - z) -
               log(1 - z) * (2polylog(3, 1 - z) + zeta(Arb(3))) + 3polylog(4, 1 - z) -
               Arb(π)^4 / 30
    end
end

function polylog_2_1_1(z::Union{Arblib.AcbOrRef,AcbSeries})
    if z isa Arblib.AcbOrRef && Arblib.contains_zero(z) && abs(z) < 1
        zᵤ = abs_ubound(Arb, z)
        return add_error(zero(z), polylog_r_div_z_bound(3, zᵤ) * zᵤ)
    else
        return Arb(π)^2 / 12 * log(1 - z)^2 +
               log(1 - z) * (polylog(3, 1 - z) + 2zeta(Arb(3))) - 3polylog(4, 1 - z) +
               Arb(π)^4 / 30
    end
end

function polylog_1_1_1_1(z::Union{Arblib.AcbOrRef,AcbSeries})
    return log(1 - z)^4 / 24
end
