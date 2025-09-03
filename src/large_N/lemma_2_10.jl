# This file contains the implementation the functions that are bounded
# in Lemma 2.10.

V_1(z) = 2polylog(1, z)

function V_1_div_z_bound(zᵤ::Arb)
    0 < zᵤ < 1 || return indeterminate(zᵤ)
    return 2polylog_r_div_z_bound(1, zᵤ)
end

V_1_log_bound_coefficients() = Arb[2]

function V_2(z)
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

function V_3_log_bound_coefficients()
    λ = λ_disc()
    return Arb[abs(λ^2 / 16 - λ + 2), abs(3λ-12)+(λ-4), 8]
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

function V_4_log_bound_coefficients()
    λ = λ_disc()
    return Arb[
        abs(λ^3/192-λ^2/8-λ/2-2)+abs(2λ*zeta(Arb(3))),
        abs(λ^2/8-2λ+4)+abs(λ^2/4-4λ+12)+abs(5λ^2/8-8λ+28),
        abs(2λ-8)+abs(6λ-24)+abs(14λ-56),
        16,
    ]
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

"""
V_log_bound_coefficients(l::Int)

Return coefficients `Cs` such that
```
abs(V(l, t)) <= sum(j -> C[j] * abs(log(1 - t))^j / factorial(j), 1:l)
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

d_0(z, t) = one(z)

d_1(z, t) = λ_disc() / 4 * log(t / z)
function d_2(z, t)
    λ = λ_disc()
    return λ^2 / 64 * log(t / z)^2 + λ / 8 * (log(t / z)^2 + 2S(2, t) - 2S(2, z))
end
function d_2_z_part(z)
    λ = λ_disc()
    return -λ / 4 * S(2, z)
end
function d_2_zt_part(z, t)
    λ = λ_disc()
    return λ^2 / 64 * log(t / z)^2 + λ / 8 * (log(t / z)^2 + 2S(2, t))
end

function d_3(z, t)
    λ = λ_disc()
    # NOTE: This is an alternative formulation that is SLIGHTLY better
    # than the version in the paper.
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + λ / 16 + 1 // 6) * log(t / z)^2 + (λ / 8 + 1) * S(2, t) -
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
function d_3_zt_part(z, t)
    λ = λ_disc()
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + λ / 16 + 1 // 6) * log(t / z)^2 + (λ / 8 + 1) * S(2, t) -
            λ / 8 * S(2, z) + S(2, conj(z))
        ) + S(3, t)
    )
end

# Part of d_3 analytic in z
function d_3_analytic_z(z, t)
    λ = λ_disc()
    return λ / 4 * (
        log(t / z) * (
            (λ^2 / 576 + (1 // 16) * λ + 1 // 6) * log(t / z)^2 +
            ((1 // 8) * λ + 1) * S(2, t) - (1 // 8) * λ * S(2, z)
        ) + S(3, t) - S(3, z)
    )
end

# Part of d_3 analytic in conj(z)
function d_3_analytic_conj_z(z, t)
    λ = λ_disc()
    return λ / 4 * log(t / z) * S(2, conj(z))
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
           λ / 4 * (S(2, t) - S(2, z)) * a * z
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
           ((λ^2 / 32 + λ / 4) * S(2, t) - (λ^2 / 32 * S(2, z) - λ / 4 * S(2, conj(z)))) *
           integral_log_z(1, z, a) +
           λ / 4 * (S(3, t) - S(3, z)) * a * z
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
    if k == 2
        d_2_z_part(z)
    elseif k == 3
        d_3_z_part(z)
    end

# Part of d only depending on z and t
d_zt_part(k::Int, z, t) =
    if k == 2
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

function integral_d_k_V_l(k::Int, l::Int, z::Acb)
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

        V_div_t * integral_d(k, z, a)
    end

    b = Arblib.contains(z, Acb(1)) ? Arb(0.999) : Arb(1)
    bz = b * z

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

        # TODO: Verify analyticity
        res_az_thin_bz_thin = if k == 0 || k == 1
            Arblib.integrate(
                az_thin,
                bz_thin,
                atol = 1e-6,
                opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
                warn_on_no_convergence = false,
            ) do t
                d(k, z, t) * V(l, t) / t
            end
        else
            # Split d into two terms, one depending only on z and one
            # depending on both z and t.

            # Part of d(k, z, t) only depending on z, so we can factor it
            # out.
            res_az_thin_bz_thin_part_1 =
                d_z_part(k, z) * Arblib.integrate(
                    az_thin,
                    bz_thin,
                    atol = 1e-6,
                    opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
                    warn_on_no_convergence = false,
                ) do t
                    V(l, t) / t
                end

            # Part of d(k, z, t) depending on both z and t.
            res_az_thin_bz_thin_part_2 = Arblib.integrate(
                az_thin,
                bz_thin,
                atol = 1e-6,
                opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
                warn_on_no_convergence = false,
            ) do t
                d_zt_part(k, z, t) * V(l, t) / t
            end

            res_az_thin_bz_thin_part_1 + res_az_thin_bz_thin_part_2
        end

        if iswide(z)
            # Add enclosures of integral from az to az_thin and from
            # z_thin to z.
            # IMPROVE: Get better enclosures of this using mean value
            # theorem?
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
            # Factor out d(k, z, t) / t from the integral. That leaves
            # ∫ V(l, t) dt from t = b * z to z.
            # The change of variables s = t / z gives
            # ∫ V(l, s * z) ds from s = b to 1.
            #
            Cs = V_log_bound_coefficients(l)
            d(k, z, t) / t * sum(1:l) do j
                Cs[j] / factorial(j) * integral_log_1mtz(Acb(1), j, b)
            end
        end
    end

    return res_0_az + res_az_bz + res_bz_z
end

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

    return t::Arblib.AcbOrRef -> begin
        M = K(convert(Acb, t))
        M.p[0] = 0
        M.p[1] = 0
        M.p[2] = 0
        M.p[3] = 0
        return (M << 4)(inv_N)
    end
end

function integral_K_4(N₀::Int, z::Acb, a::Arb)
    return zero(z) # FIXME: Implement this
end

function integral_K_4_V_1(N₀::Int, z::Acb)
    K4 = K_4(N₀, z)

    a = Arb(1e-4)
    az = a * z

    # Integrate from 0 to a * z
    res_0_az = let
        # Factor out enclosure of V(1, t) / t and integrate K_4(N₀, z, t)

        # Compute enclosure of V(1, t) / t for t in [0, a * z]
        # IMPROVE: Improving this enclosure would allow us to take a
        # larger a, which should make the integration much faster.
        azᵤ = abs_ubound(Arb, az)
        V_div_t = add_error(Acb(0), V_div_z_bound(1, azᵤ))

        V_div_t * integral_K_4(N₀, z, a)
    end

    b = Arblib.contains(z, Acb(1)) ? Arb(0.999) : Arb(1)
    bz = b * z

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

        # TODO: Verify analyticity
        # IMPROVE: Look at splitting K_4 into part depending on t and
        # part not depending on t.
        res_az_thin_bz_thin = Arblib.integrate(
            az_thin,
            bz_thin,
            atol = 1e-6,
            opts = Arblib.calc_integrate_opt_struct(0, 4_000, 0, 1, 0),
            warn_on_no_convergence = false,
        ) do t
            K4(t) * V(1, t) / t
        end

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
            # Factor out K_4(N₀, z, t) / t from the integral. That leaves
            # ∫ V(l, t) dt from t = b * z to z.
            # The change of variables s = t / z gives
            # ∫ V(l, s * z) ds from s = b to 1.
            #
            Cs = V_log_bound_coefficients(1)
            K4(t) / t * Cs[1] * integral_log_1mtz(Acb(1), 1, b)
        end
    end

    return res_0_az + res_az_bz + res_bz_z
end
