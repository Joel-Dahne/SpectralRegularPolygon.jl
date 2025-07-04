###
# Lemma 2.6
###

function lemma_2_6_c_2(; verbose = false)
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        c_2(exp(Acb(0, θ)))
    end
end

function lemma_2_6_c_3(; verbose = false)
    # FIXME: Set left bound to 0
    ArbExtras.maximum_enclosure(
        Arf(0.01),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        c_3(exp(Acb(0, θ)))
    end
end

function lemma_2_6_c_4(; verbose = false)
    # FIXME: Set left bound to 0
    ArbExtras.maximum_enclosure(
        Arf(0.01),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        maxevals = 40000,
        depth = 30,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        abs(c_4(exp(Acb(0, θ))))
    end
end

function lemma_2_6_c_5(; verbose = false)
    # FIXME: Set left bound to 0
    ArbExtras.maximum_enclosure(
        Arf(0.01),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        maxevals = 40000,
        depth = 30,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        abs(c_5(exp(Acb(0, θ))))
    end
end

function lemma_2_6_T_2(; verbose = false)
    return indeterminate(Arb)
end

function lemma_2_6_T_4(; verbose = false)
    return indeterminate(Arb)
end

function lemma_2_6_T_6(; verbose = false)
    return indeterminate(Arb)
end

function lemma_2_6(; verbose = false)
    c_2_bound = lemma_2_6_c_2(; verbose)
    c_3_bound = lemma_2_6_c_3(; verbose)
    c_4_bound = lemma_2_6_c_4(; verbose)
    c_5_bound = indeterminate(Arb) # lemma_2_6_c_5(; verbose)

    c_bounds = [c_2_bound, c_3_bound, c_4_bound, c_5_bound]

    T_2_bound = lemma_2_6_T_2(; verbose)
    T_4_bound = lemma_2_6_T_2(; verbose)
    T_6_bound = lemma_2_6_T_2(; verbose)

    T_bounds = [T_2_bound, T_4_bound, T_6_bound]

    return c_bounds, T_bounds
end

###
# Lemma 2.7
###

function lemma_2_7(; verbose = false)
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

function lemma_2_10_d_k_V_l(k::Int, l::Int; verbose = false)
    a = Arf(1e-3)

    verbose && @info "Splitting interval [0, a] and [a, π]" a

    # Compute maximum on [0, a]
    res1 = let
        z = exp(Acb(0, Arb((0, a))))
        b = exp(Acb(0, a))

        # Integrate from 0 to b
        term1 = integral_d_k_V_l_other_limit(k, l, z, b)

        # Integrate from b to z
        # TODO: Implement this
        term2 = zero(term1)

        abs(real(term1 + term2))
    end

    verbose && @info "Maximum for θ in [0, a]" res1

    # Compute maximum on [a, π]
    res2 = ArbExtras.maximum_enclosure(
        a,
        Arblib.ubound(Arb(π)),
        degree = -1,
        atol = 1e-1,
        depth_start = 4,
        abs_value = true,
        maxevals = 4096,
        threaded = true;
        verbose,
    ) do θ
        z = exp(Acb(0, θ))
        real(integral_d_k_V_l(k, l, z))
    end

    verbose && @info "Maximum for θ in [a, π]" res2

    return max(res1, res2)
end

function lemma_2_10(; verbose = false)
    kls = [
        (0, 1),
        (0, 2),
        (0, 3),
        #(0, 4),
        (1, 1),
        (1, 2),
        (1, 3),
        #(1, 4),
        (2, 1),
        (2, 2),
        (2, 3),
        (3, 1),
        (3, 2),
    ]

    res = map(kls) do (k, l)
        verbose && @info "Computing bounds for k = $k, l = $l"
        lemma_2_10_d_k_V_l(k, l; verbose)
    end

    return kls, res
end

###
# Lemma 2.13
###

function lemma_2_13(; verbose = false)
    λ = λ_disc()

    return Arb(1 / (sqrt(λ) * besselj1(sqrt(λ))), prec = 53)
end
