###
# Lemma 2.6
###

function lemma_2_6_c_2(; verbose = false)
    # Use symmetry when reflecting with real axis
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
    # Use symmetry when reflecting with real axis
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        maxevals = 10000,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        c_3(exp(Acb(0, θ)))
    end
end

function lemma_2_6_c_4(; verbose = false)
    # Use symmetry when reflecting with real axis
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        maxevals = 10000,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        abs(c_4(exp(Acb(0, θ))))
    end
end

function lemma_2_6_c_5(; verbose = false)
    # Use symmetry when reflecting with real axis
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        maxevals = 40000,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        abs(c_5(exp(Acb(0, θ))))
    end
end

function lemma_2_6_T_2(; verbose = false)
    # TODO: Implement this
end

function lemma_2_6_T_4(; verbose = false)
    # TODO: Implement this
end

function lemma_2_6_T_6(; verbose = false)
    # TODO: Implement this
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

    return res1, res2
end

###
# Lemma 2.10
###

function lemma_2_10_d_k_V_l(k::Int, l::Int; verbose = false)
    a = Arf(1e-3)

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

    verbose && @info "Maximum for θ in [0, a]" a res1
    return res1
    # Compute maximum on [a, π]
    # PROVE: That we only have to consider [0, π]
    res2 = ArbExtras.maximum_enclosure(
        a,
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        depth_start = 4,
        abs_value = true,
        maxevals = 4096,
        threaded = true;
        verbose,
    ) do θ
        z = exp(Acb(0, θ))
        real(integral_d_k_V_l(k, l, z))
    end

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
