###
# Lemma 2.6
###

function lemma_2_6_b_2(; verbose = true)
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        b_2(exp(Acb(0, θ)))
    end
end

function lemma_2_6_b_3(; verbose = true)
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        b_3(exp(Acb(0, θ)))
    end
end

function lemma_2_6_b_4(; verbose = true)
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        depth = 30,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        b_4(exp(Acb(0, θ)))
    end
end

function lemma_2_6_b_5(; verbose = true)
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        depth = 30,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        b_5(exp(Acb(0, θ)))
    end
end

function lemma_2_6_T_6(N₀::Int = 26; verbose = true)
    f = T_6_bound(N₀)

    # In practice maximum is attained at z = 1.
    maximum_estimate = abs(f(Acb(1)))
    maximum_estimate_format =
        ArbExtras.format_interval(Arblib.getinterval(maximum_estimate)...)
    verbose && @info "Enclosure at z = 1: $(maximum_estimate_format)"

    # We are happy if our enclosure is slightly larger than that.
    ubound_tol = 1.01Arblib.ubound(maximum_estimate)
    verbose && @info "Using upper tolerance: $ubound_tol"

    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        ubound_tol = ubound_tol,
        depth = 30,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        f(exp(Acb(0, θ)))
    end
end

function lemma_2_6(N₀::Int = 26; verbose = true)
    b_2_bound = lemma_2_6_b_2(; verbose)
    b_3_bound = lemma_2_6_b_3(; verbose)
    b_4_bound = lemma_2_6_b_4(; verbose)
    b_5_bound = lemma_2_6_b_5(; verbose)

    b_bounds = [b_2_bound, b_3_bound, b_4_bound, b_5_bound]

    T_6_bound = lemma_2_6_T_6(N₀; verbose)

    T_2_bound =
        b_2_bound + b_3_bound / N₀ + b_4_bound / N₀^2 + b_5_bound / N₀^3 + T_6_bound / N₀^4
    T_4_bound = b_4_bound + b_5_bound / N₀ + T_6_bound / N₀^2

    T_bounds = [T_2_bound, T_4_bound, T_6_bound]

    return b_bounds, T_bounds
end

###
# Lemma 2.7
###

function lemma_2_7(; verbose = true)
    b = Arf(5)

    res1 = ArbExtras.maximum_enclosure(
        g_dw3,
        Arf(0),
        b,
        rtol = 1e-3,
        abs_value = true,
        threaded = true;
        verbose,
    )

    # Bound on [b, Inf]
    res2 = Arb("0.7858") * λ_disc()^(3 // 2 - 1 // 6) * Arb(b)^(-1 // 3)

    verbose && @info "Bound on [$b, Inf] is $res2"

    return max(res1, res2)
end

###
# Lemma 2.10
###

function lemma_2_10_d_k_V_l(k::Int, l::Int; ubound_tol::Arb = Arb(0), verbose = true)
    a = Arf(1e-3)

    verbose && @info "Splitting interval [0, a] and [a, π]" a

    # Compute maximum on [0, a]
    res_0_a = Arb(0)
    # TODO: Implement this
    # res_0_a = let
    #     z = exp(Acb(0, Arb((0, a))))
    #     b = exp(Acb(0, a))

    #     # Integrate from 0 to b
    #     term_0_b = integral_d_k_V_l_other_limit(k, l, z, b)

    #     # Integrate from b to z
    #     # TODO: Implement this
    #     term_b_z = zero(term_0_b)

    #     abs(real(term_0_b + term_b_z))
    # end

    verbose && @info "Maximum for θ in [0, a]" res_0_a

    verbose && @info "Computing maximum on [a, π]" ubound_tol

    # Compute maximum on [a, π]
    res_a_π = ArbExtras.maximum_enclosure(
        a,
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        depth_start = 4,
        abs_value = true,
        depth = 30,
        maxevals = 2000,
        threaded = true;
        ubound_tol,
        verbose,
    ) do θ
        z = exp(Acb(0, θ))
        real(integral_d_k_V_l(k, l, z))
    end

    verbose && @info "Maximum for θ in [a, π]" res_a_π

    return max(res_0_a, res_a_π)
end

function lemma_2_10(; verbose = true)
    kls = [
        #(0, 1),
        #(0, 2),
        #(0, 3),
        #(0, 4),
        #(1, 1),
        #(1, 2),
        #(1, 3),
        (1, 4), # TODO
        #(2, 1),
        #(2, 2),
        (2, 3), # TODO
        #(3, 1),
        (3, 2), # TODO
    ]

    res = map(kls) do (k, l)
        verbose && @info "Computing bounds for k = $k, l = $l"
        lemma_2_10_d_k_V_l(k, l; verbose)
    end

    # TODO: Bound integral with K_4

    return kls, res
end

###
# Lemma 2.13
###

function lemma_2_13(; verbose = true)
    λ = λ_disc()

    return Arb(1 / (sqrt(λ) * besselj1(sqrt(λ))), prec = 53)
end
