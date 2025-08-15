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

function lemma_2_16(; verbose = true)

end

###
# Lemma 2.19
###

function lemma_2_19(; verbose = true)

end

###
# Lemma 2.20
###

function lemma_2_20(; verbose = true)

end
