"""
    rgamma(x)

Compute the reciprocal gamma function, defined by `rgamma(x) = 1 /
gamma(x)`.
"""
rgamma(x::Union{Arb,ArbSeries,AcbSeries}) = Arblib.rgamma!(zero(x), x)

hypgeom0f1_regularized(a::Arb, z::Arb) = Arblib.hypgeom_0f1!(zero(z), a, z, 1)
hypgeom0f1_regularized(a::Acb, z::Acb) = Arblib.hypgeom_0f1!(zero(z), a, z, 1)

function hypgeom0f1_regularized(a::Arb, z::ArbSeries)
    z0 = z[0]
    res = zero(z)

    for i = 0:Arblib.degree(z)
        res[i] = hypgeom0f1_regularized(a + i, z0) / factorial(BigInt(i))
    end

    return ArbExtras.compose_zero!(res, res, z)
end

function hypgeom0f1_regularized(a::Acb, z::AcbSeries)
    z0 = z[0]
    res = zero(z)

    for i = 0:Arblib.degree(z)
        res[i] = hypgeom0f1_regularized(a + i, z0) / factorial(BigInt(i))
    end

    return ArbExtras.compose_zero!(res, res, z)
end

hypgeom2f1(a::Acb, b::Acb, c::Acb, z::Acb) =
    Arblib.hypgeom_2f1!(zero(z), a, b, c, z, flags = 0)
