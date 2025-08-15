###
# Lemma 2.19
###

function lemma_2_19(N₀::Int; verbose = true)
    inv_N = Arb((0, 1 // N₀))

    lhs = λ_approx(inv_N) / (1 - ϵ_prime(inv_N))
    rhs = λ₂_disc() / Arb("1.01")^2

    return lhs, rhs
end

###
# Lemma 2.20
###

function lemma_2_20_part1(N₀::Int; verbose = true)
    # Check that first three derivatives are zero at zero
    expansion_zero = λ_sup_m_λ_inf(ArbSeries((0, 1), degree = 3))
    @assert iszero(expansion_zero[0])
    @assert iszero(expansion_zero[1])
    @assert iszero(expansion_zero[2])
    @warn "TODO: Implement check for third term" iszero(expansion_zero[3])

    # Check that fourth derivative is positive
    return ArbExtras.minimum_enclosure(
        ArbExtras.derivative_function(λ_sup_m_λ_inf, 4),
        Arf(0),
        ubound(Arb(1 // N₀));
        verbose,
    )
end

function lemma_2_20_part2(N₀::Int; verbose = true)
    # Check that first three derivatives are zero at zero
    expansion_zero = q_sup_m_q_inf(ArbSeries((0, 1), degree = 4))
    @warn "TODO: Implement check for first term" iszero(expansion_zero[3])
    @assert iszero(expansion_zero[1])
    @assert iszero(expansion_zero[2])
    @warn "TODO: Implement check for third term" iszero(expansion_zero[3])
    @warn "TODO: Implement check for fourth term" iszero(expansion_zero[4])

    a = Arf(1e-3)

    # Check that fourth derivative is positive on [0, a]
    derivative_0_a = ArbExtras.minimum_enclosure(
        ArbExtras.derivative_function(q_sup_m_q_inf, 5),
        Arf(0),
        a;
        verbose,
    )

    res_a_inv_N =
        ArbExtras.minimum_enclosure(q_sup_m_q_inf, a, ubound(Arb(1 // N₀)); verbose)

    return derivative_0_a, res_a_inv_N
end
