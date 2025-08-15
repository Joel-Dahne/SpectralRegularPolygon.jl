###
# Lemma 2.15
###

function lemma_2_15(N₀; verbose = true)
    A = Arb(π)
    inv_N = Arb((0, 1 // N₀))

    # TODO: Prove monotonicity?
    I = sqrt(A / (π * sinc(inv_N) / cospi(inv_N)))

    # TODO: Prove monotonicity?
    C = sqrt(A / (π * sinc(2inv_N)))

    I > Arb("0.99") || throw(ErrorException("bound I > 0.99 doesn't hold"))

    C < Arb("1.01") || throw(ErrorException("bound C < 1.01 doesn't hold"))

    return I, C
end

###
# Lemma 2.16
###

function lemma_2_16_V_l(l::Int, N₀::Int; verbose = true)
    # FIXME: Add correct lower bound for abs(z) here
    R = Arb("0.99")

    inv_N = Arb((0, 1 // N₀))
    z_pow_N = add_error(Acb(0), R^N₀)

    f_prime = c_N(N₀) / (1 - z_pow_N)^2inv_N

    V(l, z_pow_N) / f_prime
end

function lemma_2_16(N₀::Int; verbose = true)
    return lemma_2_16_V_l.(1:4, N₀; verbose)
end
