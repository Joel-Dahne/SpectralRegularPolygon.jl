hypgeom0f1_regularized(a::Arb, z::Arb) = Arblib.hypgeom_0f1!(zero(z), a, z, 1)

function hypgeom0f1_regularized(a::Arb, z::ArbSeries)
    z0 = z[0]
    res = zero(z)

    for i = 0:Arblib.degree(z)
        res[i] = hypgeom0f1_regularized(a + i, z0) / factorial(Arb(i))
    end

    return ArbExtras.compose_zero!(res, res, z)
end
