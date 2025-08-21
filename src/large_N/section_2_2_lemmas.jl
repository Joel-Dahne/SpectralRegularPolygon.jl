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

function lemma_2_16_V_l(l::Int, N₀::Int, R₂::Arb; verbose = true)
    inv_N = Arb((0, 1 // N₀))
    z_pow_N = add_error(Acb(0), R₂^N₀)

    f_prime = c_N(N₀) / (1 - z_pow_N)^2inv_N

    abs(V(l, z_pow_N) / f_prime)
end

function lemma_2_16(N₀::Int; verbose = true)
    R = Arb("0.95")
    R₂ = Arb("0.951")

    N_max = 128

    # Verify that R₂ works as a lower bound
    res1 = map(N₀:N_max) do N
        let inv_N = Acb(1 // N)
            _c_N(real(inv_N)) * R₂ * real(hypgeom2f1(2inv_N, inv_N, 1 + inv_N, Acb(-R₂^N)))
        end
    end

    res2 = let inv_N = Acb(Arb((0, 1 // N_max)))
        c_N(N₀) * R₂ * real(hypgeom2f1(2inv_N, inv_N, 1 + inv_N, Acb(Arb((-R₂^N_max, 0)))))
    end

    all(R .< res1) && (R < res2) || throw(ErrorException("verification of R₂ failed"))

    return lemma_2_16_V_l.(1:4, N₀, R₂; verbose)
end
