# This file contains the implementation the functions that are bounded
# in Lemma 2.10.

function V_1(z; analytic::Bool = false)
    # Check analyticity when requested. It has a branch cut for z in
    # [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(z)) &&
       Arblib.contains_nonnegative(Arblib.realref(z) - 1)

        return indeterminate(z)
    end

    2polylog(1, z)
end

function V_1_div_z_bound(zᵤ::Arb)
    0 < zᵤ < 1 || return indeterminate(zᵤ)
    return 2polylog_r_div_z_bound(1, zᵤ)
end

V_1_log_bound_coefficients() = Arb[2]

function V_2(z; analytic::Bool = false)
    # Check analyticity when requested. It has a branch cut for z in
    # [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(z)) &&
       Arblib.contains_nonnegative(Arblib.realref(z) - 1)

        return indeterminate(z)
    end

    λ = λ_disc()
    return (λ / 2 - 2) * polylog(2, z) + 4polylog_1_1(z)
end

function V_2_div_z_bound(zᵤ::Arb)
    0 < zᵤ < 1 || return indeterminate(zᵤ)
    λ = λ_disc()
    return abs(λ / 2 - 2) * polylog_r_div_z_bound(1, zᵤ) + 4polylog_r_div_z_bound(2, zᵤ)
end

function V_2_log_bound_coefficients()
    λ = λ_disc()
    return Arb[abs(λ / 2 - 2), 4]
end

function V_3(z; analytic::Bool = false)
    # Check analyticity when requested. It has a branch cut for z in
    # [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(z)) &&
       Arblib.contains_nonnegative(Arblib.realref(z) - 1)

        return indeterminate(z)
    end

    λ = λ_disc()
    return (λ^2 / 16 - λ + 2) * polylog(3, z) +
           (3λ - 12) * polylog_1_2(z) +
           (λ - 4) * polylog_2_1(z) +
           8polylog_1_1_1(z)
end

function V_3_div_z_bound(zᵤ::Arb)
    0 < zᵤ < 1 || return indeterminate(zᵤ)
    λ = λ_disc()
    return abs(λ^2 / 16 - λ + 2) * polylog_r_div_z_bound(1, zᵤ) +
           (abs(3λ - 12) + abs(λ - 4)) * polylog_r_div_z_bound(2, zᵤ) +
           8polylog_r_div_z_bound(3, zᵤ)
end

function V_3_log_bound_coefficients()
    λ = λ_disc()
    return Arb[abs(λ^2 / 16 - λ + 2), abs(3λ-12)+abs(λ-4), 8]
end

function V_4(z; analytic::Bool = false)
    # Check analyticity when requested. It has a branch cut for z in
    # [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(z)) &&
       Arblib.contains_nonnegative(Arblib.realref(z) - 1)

        return indeterminate(z)
    end

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
    return (abs(λ^3 / 192 - λ^2 / 8 - λ / 2 - 2) + abs(2λ * zeta(Arb(3)))) *
           polylog_r_div_z_bound(1, zᵤ) +
           (abs(λ^2 / 8 - 2λ + 4) + abs(λ^2 / 4 - 4λ + 12) + abs(5λ^2 / 8 - 8λ + 28)) *
           polylog_r_div_z_bound(2, zᵤ) +
           (abs(2λ - 8) + abs(6λ - 24) + abs(14λ - 56)) * polylog_r_div_z_bound(3, zᵤ) +
           16polylog_r_div_z_bound(4, zᵤ)
end

function V_4_log_bound_coefficients()
    λ = λ_disc()
    return Arb[
        abs(λ^3/192-λ^2/8-λ/2-2)+abs(2λ*zeta(Arb(3))),
        abs(λ^2/8-2λ+4)+abs(λ^2/4-4λ+12)+abs(5λ^2/8-8λ+28),
        abs(2λ-8)+abs(6λ-24)+abs(14λ-56),
        16,
    ]
end

V(l::Int, z::Arblib.AcbOrRef; analytic::Bool = false) =
    if l == 1
        V_1(z; analytic)
    elseif l == 2
        V_2(z; analytic)
    elseif l == 3
        V_3(z; analytic)
    elseif l == 4
        V_4(z; analytic)
    end

"""
    V_div_z_bound(l::Int, zᵤ::Arb)

Return a constant `D` such that
```
abs(V(l, z)) <= D * z
```
for all `|z| < zᵤ < 1`.

This is based on Lemma A.2 in the paper.
"""
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

"""
V_log_bound_coefficients(l::Int)

Return coefficients `Cs` such that
```
abs(V(l, t)) <= sum(j -> C[j] * abs(log(1 - t))^j / factorial(j), 1:l)
```
FIXME: Rewrite this to instead give the bound
```
abs(V(l, t)) <= sum(j -> C[j] * log(1 - abs(t))^j / factorial(j), 1:l)
```
"""
V_log_bound_coefficients(l::Int) =
    if l == 1
        V_1_log_bound_coefficients()
    elseif l == 2
        V_2_log_bound_coefficients()
    elseif l == 3
        V_3_log_bound_coefficients()
    elseif l == 4
        V_4_log_bound_coefficients()
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

# Part of d only depending on z
d_z_part(k::Int, z) =
    if k == 2
        d_2_z_part(z)
    elseif k == 3
        d_3_z_part(z)
    end

# Part of d only depending on z and t
d_zt_part(k::Int, z, t; analytic::Bool = false) =
    if k == 2
        d_2_zt_part(z, t; analytic)
    elseif k == 3
        d_3_zt_part(z, t; analytic)
    end

"""
integral_d(k::Int, z::Acb, a::Arb)

Compute a bound for the integral of
```
∫ abs(d(l, z, t)) abs(dt)
```
taken from `0` to `a * z`.
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
        # IMPROVE: Improving this enclosure would allow us to take a
        # larger a, which should make the integration much faster.
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
            Cs = V_log_bound_coefficients(l)
            integral_V_bound = sum(1:l) do j
                Cs[j] / factorial(j) * abs(integral_log_1mtz(Acb(1), j, b))
            end

            add_error(Acb(0), d_k_div_t_bound * integral_V_bound)
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

TODO: Write about the above in the paper.
TODO: Should we use that `abs(z) = 1`?
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

    return t::Acb -> begin
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
    K = K_model(N₀, z)

    return (t::Arblib.AcbOrRef; analytic::Bool = false) -> begin
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
    K4 = K_4(N₀, z)

    a = Arb(1e-4)
    b = Arblib.contains(z, Acb(1)) ? Arb(1) - 1e-6 : Arb(1)
    az = a * z
    bz = b * z

    # Integrate from 0 to a * z
    res_0_az = let
        # V(l, t) / t is bounded near t = 0, hence we can factor out a
        # bound of it and integrate only K_4(z, t).


        # Compute bound of V(l, t) / t for t in [0, a * z]
        # IMPROVE: Improving this enclosure would allow us to take a
        # larger a, which should make the integration much faster.
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
            # IMPROVE: Get better enclosures of this using mean value
            # theorem?
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
            Cs = V_log_bound_coefficients(1)
            # This sum only has one term, but we write it like this to
            # make the form clearer.
            integral_V_bound = sum(1:1) do j
                Cs[j] / factorial(j) * abs(integral_log_1mtz(Acb(1), j, b))
            end

            add_error(Acb(0), K_4_div_t_bound * integral_V_bound)
        end
    end

    return res_0_az + res_az_bz + res_bz_z
end
