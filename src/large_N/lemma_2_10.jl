# This file contains the implementation the functions that are bounded
# in Lemma 2.10.

"""
    d(k::Int, z::Arblib.AcbOrRef, t::Arblib.AcbOrRef; analytic::Bool = false)

Compute an enclosure of `d_k(z, t)` from Equation REF(21) in the paper.

If `analytic` is true, then return an indeterminate value if `t`
overlaps a branch cut of the function.
"""
d(k::Int, z::Arblib.AcbOrRef, t::Arblib.AcbOrRef; analytic::Bool = false) =
    if k == 0
        d_0(z, t; analytic)
    elseif k == 1
        d_1(z, t; analytic)
    elseif k == 2
        d_2(z, t; analytic)
    elseif k == 3
        d_3(z, t; analytic)
    end

"""
    d_z_part(k::Int, z::Arblib.AcbOrRef)

Compute an enclosure of the terms of [`d`](@ref) which only depend on
`z`.
"""
d_z_part(k::Int, z::Arblib.AcbOrRef) =
    if k == 2
        d_2_z_part(z)
    elseif k == 3
        d_3_z_part(z)
    end

"""
    d_zt_part(k::Int, z::Arblib.AcbOrRef, t::Arblib.AcbOrRef; analytic::Bool = false)

Compute an enclosure of the terms of [`d`](@ref) which depend on both
`z` and `t`.

If `analytic` is true, then return an indeterminate value if `t`
overlaps a branch cut of the function.
"""
d_zt_part(k::Int, z::Arblib.AcbOrRef, t::Arblib.AcbOrRef; analytic::Bool = false) =
    if k == 2
        d_2_zt_part(z, t; analytic)
    elseif k == 3
        d_3_zt_part(z, t; analytic)
    end

d_0(z, t; analytic::Bool = false) = one(z)

function d_1(z, t; analytic::Bool = false)
    # Compute log(t / z), checking analyticity when requested
    log_t_div_z = Arblib.log_analytic!(zero(t), t / z, analytic)
    isfinite(log_t_div_z) || return indeterminate(t)

    λ_disc() / 4 * log_t_div_z
end

function d_2(z, t; analytic::Bool = false)
    # Check analyticity of S(2, t) when requested. It has a branch cut
    # for t in [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(t)) &&
       Arblib.contains_nonnegative(Arblib.realref(t) - 1)

        return indeterminate(t)
    end

    λ = λ_disc()
    # Compute log(t / z), checking analyticity when requested
    log_t_div_z = Arblib.log_analytic!(zero(t), t / z, analytic)
    isfinite(log_t_div_z) || return indeterminate(t)

    return λ^2 / 64 * log_t_div_z^2 + λ / 8 * (log_t_div_z^2 + 2S(2, t) - 2S(2, z))
end
function d_2_z_part(z)
    λ = λ_disc()
    return -λ / 4 * S(2, z)
end
function d_2_zt_part(z, t; analytic::Bool = false)
    # Check analyticity of S(2, t) when requested. It has a branch cut
    # for t in [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(t)) &&
       Arblib.contains_nonnegative(Arblib.realref(t) - 1)

        return indeterminate(t)
    end

    λ = λ_disc()
    # Compute log(t / z), checking analyticity when requested
    log_t_div_z = Arblib.log_analytic!(zero(t), t / z, analytic)
    isfinite(log_t_div_z) || return indeterminate(t)

    return λ^2 / 64 * log_t_div_z^2 + λ / 8 * (log_t_div_z^2 + 2S(2, t))
end

function d_3(z, t; analytic::Bool = false)
    # Check analyticity of S(2, t) and S(3, t) when requested. It has
    # a branch cut for t in [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(t)) &&
       Arblib.contains_nonnegative(Arblib.realref(t) - 1)

        return indeterminate(t)
    end

    λ = λ_disc()
    # Compute log(t / z), checking analyticity when requested
    log_t_div_z = Arblib.log_analytic!(zero(t), t / z, analytic)
    isfinite(log_t_div_z) || return indeterminate(t)

    # NOTE: This is an alternative formulation that is SLIGHTLY better
    # than the version in the paper.
    return λ / 4 * (
        log_t_div_z * (
            (λ^2 / 576 + λ / 16 + 1 // 6) * log_t_div_z^2 + (λ / 8 + 1) * S(2, t) -
            λ / 8 * S(2, z) + S(2, conj(z))
        ) + S(3, t) - S(3, z)
    )
    # This is the version in the paper
    #return λ^3 / 2304 * log(t / z)^3 +
    #       λ^2 / 64 * log(t / z) * (log(t / z)^2 + 2S(2, t) - 2S(2, z)) +
    #       λ / 24 * (
    #           log(t / z)^3 +
    #           6log(t / z) * S(2, t) +
    #           6log(t / z) * S(2, conj(z)) +
    #           6S(3, t) - 6S(3, z)
    #       )
end
function d_3_z_part(z)
    λ = λ_disc()
    return -λ / 4 * S(3, z)
end
function d_3_zt_part(z, t; analytic::Bool = false)
    # Check analyticity of S(2, t) and S(3, t) when requested. It has
    # a branch cut for t in [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(t)) &&
       Arblib.contains_nonnegative(Arblib.realref(t) - 1)

        return indeterminate(t)
    end

    λ = λ_disc()
    # Compute log(t / z), checking analyticity when requested
    log_t_div_z = Arblib.log_analytic!(zero(t), t / z, analytic)
    isfinite(log_t_div_z) || return indeterminate(t)

    return λ / 4 * (
        log_t_div_z * (
            (λ^2 / 576 + λ / 16 + 1 // 6) * log_t_div_z^2 + (λ / 8 + 1) * S(2, t) -
            λ / 8 * S(2, z) + S(2, conj(z))
        ) + S(3, t)
    )
end

"""
    integral_d(k::Int, z::Acb, a::Arb)

Compute a bound for the integral of
```
∫ abs(d(l, z, t)) abs(dt)
```
taken from `0` to `a * z`.

This is based on Lemma REF(A.3) in the paper.
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

function integral_d_0(z::Acb, a::Arb)
    return a * abs(z)
end

function integral_d_1(z::Acb, a::Arb)
    λ = λ_disc()
    return λ / 4 * abs(z) * abs(integral_log(1, a))
end

function integral_d_2(z::Acb, a::Arb)
    λ = λ_disc()
    t = Arblib.union(zero(z), a * z)
    return (λ^2 / 64 + λ / 8) * abs(z) * abs(integral_log(2, a)) +
           λ / 4 * abs(S(2, t) - S(2, z)) * a * abs(z)
end

function integral_d_3(z::Acb, a::Arb)
    λ = λ_disc()
    t = Arblib.union(zero(z), a * z)
    return (λ^3 / 2304 + λ^2 / 64 + λ / 24) * abs(z) * abs(integral_log(3, a)) +
           (λ^2 / 32 * abs(S(2, t) - S(2, z)) + λ / 4 * abs(S(2, t) + S(2, conj(z)))) *
           abs(z) *
           integral_log(1, a) +
           λ / 4 * abs(S(3, t) - S(3, z)) * a * abs(z)
end

"""
    integral_d_k_V_l(k::Int, l::Int, z::Acb)

Compute an enclosure of the integral
```
∫ d_k(z, t) * V_l(t) dt
```
from `0` to `z`, which is bounded in Lemma REF(2.10) in the paper.
"""
function integral_d_k_V_l(k::Int, l::Int, z::Acb)
    a = Arb(1e-4)
    b = Arblib.contains(z, Acb(1)) ? Arb(0.999) : Arb(1)
    az = a * z
    bz = b * z

    # Integrate from 0 to a * z
    res_0_az = let
        # V(l, t) / t is bounded near t = 0, hence we can factor out a
        # bound of it and integrate only d(k, z, t).

        # Compute bound of V(l, t) / t for t in [0, a * z]
        azᵤ = abs_ubound(Arb, az)
        D_l_a = V_div_z_bound(l, azᵤ)

        add_error(Acb(0), D_l_a * abs(integral_d(k, z, a)))
    end

    # Integrate from a * z to b * z
    res_az_bz = let
        # Arblib.integrate doesn't handle wide integration limits very
        # well. For that reason we integrate from the midpoint of the
        # endpoints if z is wide, and add the remaining part later.
        az_thin, bz_thin = if iswide(z)
            midpoint(Acb, az), midpoint(Acb, bz)
        else
            az, bz
        end

        res_az_thin_bz_thin = if k == 0 || k == 1
            Arblib.integrate(
                (t; analytic) -> d(k, z, t; analytic) * V(l, t; analytic) / t,
                az_thin,
                bz_thin,
                check_analytic = true,
                atol = 1e-6,
                opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
                warn_on_no_convergence = false,
            )
        else
            # Split d into two terms, one depending only on z and one
            # depending on both z and t.

            # Part of d(k, z, t) only depending on z, so we can factor it
            # out.
            res_az_thin_bz_thin_part_1 =
                d_z_part(k, z) * Arblib.integrate(
                    (t; analytic) -> V(l, t; analytic) / t,
                    az_thin,
                    bz_thin,
                    check_analytic = true,
                    atol = 1e-6,
                    opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
                    warn_on_no_convergence = false,
                )

            # Part of d(k, z, t) depending on both z and t.
            res_az_thin_bz_thin_part_2 = Arblib.integrate(
                (t; analytic) -> d_zt_part(k, z, t; analytic) * V(l, t; analytic) / t,
                az_thin,
                bz_thin,
                check_analytic = true,
                atol = 1e-6,
                opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
                warn_on_no_convergence = false,
            )

            res_az_thin_bz_thin_part_1 + res_az_thin_bz_thin_part_2
        end

        if iswide(z)
            # Add enclosures of integral from az to az_thin and from
            # z_thin to z.
            (res_az_thin_bz_thin) +
            (az_thin - az) * d(k, z, az) * V(l, az) / az +
            (bz - bz_thin) * d(k, z, bz) * V(l, bz) / bz
        else
            res_az_thin_bz_thin # We integrated everything
        end
    end

    res_bz_z = if isone(b)
        zero(res_az_bz) # We integrated everything
    else
        let t = Arblib.union(bz, z)
            # d(k, z, t) / t is bounded on the interval of
            # integration, we factor out a bound for it.
            d_k_div_t_bound = abs_ubound(d(k, z, t) / t)

            # To bound ∫ V(l, t) dt from t = b * z to z we use
            # V_log_bound_coefficients to get a bound in terms of
            # logarithms, that are then explicitly integrated.
            C = V_log_bound_coefficients(l)
            integral_V_l_bound = sum(1:l) do j
                C[j] / factorial(j) * abs(integral_log_1mtz(Acb(1), j, b))
            end

            add_error(Acb(0), d_k_div_t_bound * integral_V_l_bound)
        end
    end

    return res_0_az + res_az_bz + res_bz_z
end

"""
    K_model(N₀::Int, z::Acb)

Compute an [`ArbTaylorModel`](@ref) of [`K`](@ref) in `inv(N)` that
is valid for all `N >= N₀`. It is computed with a remainder term of
degree 5.

The paper uses the formula
```
K(z, t) = besselj0(sqrt(ρ) * abs(z)^(1 / N) * sqrt(F_N(conj(z)) * (F_N(z) - (t / z)^(1 / N) * F_N(t))))
```
To avoid the square root we use that
```
besselj0(x) = hypgeom0f1_regularized(1, -x^2 / 4)
```
This gives us the formula
```
K(z, t) = hypgeom0f1_regularized(1, -ρ * abs(z)^(2 / N) * F_N(conj(z)) * (F_N(z) - (t / z)^(1 / N) * F_N(t)) / 4)
```
"""
function K_model(N₀::Int, z::Acb)
    inv_N = Arb((0, 1 // N₀))

    # It is really important to get a good enclosure of ρ, so we
    # compute a higher order Taylor model and then truncate it.
    ρ_model = truncate(
        ArbTaylorModel(inv_N -> c_N(inv_N)^2 * λ_app(inv_N), inv_N, Arb(0), degree = 10),
        degree = 5,
    )

    abs_z_pow_2inv_N_model = ArbTaylorModel(inv_N, Arb(0), degree = 5) do inv_N
        abs(z)^2inv_N
    end

    F_N_z_model = F_N_model(N₀, z)
    F_N_conj_z_model = F_N_model(N₀, conj(z))

    finite =
        isfinite(ρ_model) &&
        isfinite(abs_z_pow_2inv_N_model) &&
        isfinite(F_N_z_model) &&
        isfinite(F_N_conj_z_model)

    return finite,
    t::Acb -> begin
        t_div_z_pow_inv_N_model = AcbTaylorModel(inv_N, Arb(0), degree = 5) do inv_N
            (t / z)^inv_N
        end

        F_N_t_model = F_N_model(N₀, t)

        return compose(
            x -> hypgeom0f1_regularized(Acb(1), -x / 4),
            ρ_model *
            abs_z_pow_2inv_N_model *
            F_N_conj_z_model *
            (F_N_z_model - t_div_z_pow_inv_N_model * F_N_t_model),
        )
    end
end

function K_4(N₀::Int, z::Acb)
    inv_N = Arb((0, 1 // N₀))
    finite, K = K_model(N₀, z)

    return finite,
    (t::Arblib.AcbOrRef; analytic::Bool = false) -> begin
        if analytic
            # Check if t overlaps the branch cut. The branch cut is
            # for either t / z lying on the negative real axis or for
            # t in [1, Inf].

            # Check if t lies in [1, Inf]
            if Arblib.contains_zero(Arblib.imagref(t)) &&
               Arblib.contains_nonnegative(Arblib.realref(t) - 1)

                return indeterminate(t)
            end

            # Check if t / z lies on negative real axis
            t_div_z = t / z
            if Arblib.contains_zero(Arblib.imagref(t_div_z)) &&
               Arblib.contains_nonpositive(Arblib.realref(t_div_z))

                return indeterminate(t)
            end
        end

        M = K(convert(Acb, t))

        # Truncate it to degree 3 and take remainder term
        return remainder(truncate(M, degree = 3))
    end
end

"""
    integral_K_4(N₀::Int, z::Acb, a::Arb)

Compute an upper bound of the absolute value of the integral of `K_4`
from `0` to `a * z`. It uses the bound
```
abs(K_4(z, t)) <=
    L_0
    + L_1 * abs(log(t / z))
    + L_2 * abs(log(t / z))^2
    + L_3 * abs(log(t / z))^3
    + L_4 * abs(log(t / z))^4
    + L_6 * abs(log(t / z))^6
```
and integrates termwise.
"""
function integral_K_4_bound(N₀::Int, z::Acb, a::Arb)
    return Arb(C_L_0) * abs(integral_log(0, a)) +
           Arb(C_L_1) * abs(integral_log(1, a)) +
           Arb(C_L_2) * abs(integral_log(2, a)) +
           Arb(C_L_3) * abs(integral_log(3, a)) +
           Arb(C_L_4) * abs(integral_log(4, a)) +
           Arb(C_L_6) * abs(integral_log(6, a))
end

function integral_K_4_V_1(N₀::Int, z::Acb)
    finite_K4, K4 = K_4(N₀, z)

    finite_K4 || return indeterminate(z)

    a = Arb(1e-4)
    b = Arblib.contains(z, Acb(1)) ? Arb(1) - 1e-6 : Arb(1)
    az = a * z
    bz = b * z

    # Integrate from 0 to a * z
    res_0_az = let
        # V(l, t) / t is bounded near t = 0, hence we can factor out a
        # bound of it and integrate only K_4(z, t).


        # Compute bound of V(l, t) / t for t in [0, a * z]
        azᵤ = abs_ubound(Arb, az)
        D_1_a = V_div_z_bound(1, azᵤ)

        add_error(Acb(0), D_1_a * integral_K_4_bound(N₀, z, a))
    end

    # Integrate from a * z to b * z
    res_az_bz = let
        # Arblib.integrate doesn't handle wide integration limits very
        # well. For that reason we integrate from the midpoint of the
        # endpoints if z is wide, and add the remaining part later.
        az_thin, bz_thin = if iswide(z)
            midpoint(Acb, az), midpoint(Acb, bz)
        else
            az, bz
        end

        res_az_thin_bz_thin = Arblib.integrate(
            (t; analytic) -> K4(t; analytic) * V(1, t; analytic) / t,
            az_thin,
            bz_thin,
            check_analytic = true,
            atol = 1e-6,
            opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
            warn_on_no_convergence = false,
        )

        if iswide(z)
            # Add enclosures of integral from az to az_thin and from
            # z_thin to z.
            (res_az_thin_bz_thin) +
            (az_thin - az) * K4(az) * V(1, az) / az +
            (bz - bz_thin) * K4(bz) * V(1, bz) / bz
        else
            res_az_thin_bz_thin # We integrated everything
        end
    end

    res_bz_z = if isone(b)
        zero(res_az_bz) # We integrated everything
    else
        let t = Arblib.union(bz, z)
            # K_4(z, t) / t is bounded on the interval of
            # integration, we factor out a bound for it.
            K_4_div_t_bound = abs_ubound(K4(t) / t)

            # To bound ∫ V(l, t) dt from t = b * z to z we use
            # V_log_bound_coefficients to get a bound in terms of
            # logarithms, that are then explicitly integrated.
            C = V_log_bound_coefficients(l)
            integral_V_l_bound = sum(1:l) do j
                C[j] / factorial(j) * abs(integral_log_1mtz(Acb(1), j, b))
            end

            add_error(Acb(0), K_4_div_t_bound * integral_V_l_bound)
        end
    end

    return res_0_az + res_az_bz + res_bz_z
end
