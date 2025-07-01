"""
    iswide(x; cutoff = 10)

Return true if `x` is wide in the meaning that the effective relative
accuracy of `x` measured in bits is more than `cutoff` lower than it's
precision. For `x` not of type `Arb` or `Acb` this always return
`false`. For `x` of type `ArbSeries` or `AcbSeries` it checks the
first coefficient.
"""
iswide(x::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}; cutoff = 10) =
    Arblib.rel_accuracy_bits(x) < precision(x) - cutoff
iswide(x::Union{ArbSeries,AcbSeries}; cutoff = 10) = iswide(Arblib.ref(x, 0); cutoff)
iswide(::Number; cutoff = 10) = false

function mean_value_theorem_bound(f, z::Arblib.AcbOrRef)
    if iswide(z)
        df_z = f(AcbSeries((z, 1)))[1]

        if isfinite(df_z)
            z0 = Arblib.midpoint(Acb, z)
            return Arblib.add_error(f(z0), abs(df_z) * abs(z - z0))
        else
            return f(z)
        end
    else
        return f(z)
    end
end
