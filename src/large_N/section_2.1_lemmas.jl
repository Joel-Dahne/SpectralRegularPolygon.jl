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
        abs(c_2(exp(Acb(0, θ))))
    end
end

function lemma_2_6_c_3(; verbose = false)
    # Use symmetry when reflecting with real axis
    ArbExtras.maximum_enclosure(
        Arf(0),
        Arblib.ubound(Arb(π)),
        degree = -1,
        rtol = 1e-3,
        maxevals = 100000,
        abs_value = true,
        threaded = true;
        verbose,
    ) do θ
        abs(c_3(exp(Acb(0, θ))))
    end
end

function lemma_2_6_c_4(; verbose = false)
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
        abs(c_3(exp(Acb(0, θ))))
    end
end

function lemma_2_6_c_5(; verbose = false)
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
        abs(c_3(exp(Acb(0, θ))))
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

###
# Lemma 2.13
###

function lemma_2_13(; verbose = false)
    λ = λ_disc()

    return Arb(1 / (sqrt(λ) * besselj1(sqrt(λ))), prec = 53)
end
