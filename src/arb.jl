# There is no version of this in Arblib
function Arblib.indeterminate!(x::Union{ArbSeries,AcbSeries})
    for i = 0:Arblib.degree(x)
        Arblib.indeterminate!(Arblib.ref(x, i))
    end
    # Since we manually set the coefficients of the polynomial we
    # need to also manually set the degree.
    Arblib.cstruct(x).length = Arblib.degree(x) + 1
    return x
end

"""
    indeterminate(x)

Construct an indeterminate version of `x`.
"""
indeterminate(x::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}) = Arblib.indeterminate!(zero(x))
indeterminate(::Type{T}) where {T<:Union{Arb,Acb}} = Arblib.indeterminate!(zero(T))
indeterminate(x::Union{ArbSeries,AcbSeries}) = Arblib.indeterminate!(zero(x))
indeterminate(::Type{T}) where {T<:AbstractFloat} = convert(T, NaN)
indeterminate(::Type{Complex{T}}) where {T<:AbstractFloat} =
    convert(Complex{T}, complex(NaN, NaN))
indeterminate(x) = indeterminate(typeof(x))

function Base.hypot(x::ArbSeries, y::ArbSeries)
    res = x^2
    Arblib.add_series!(res, res, y^2, length(res))
    Arblib.sqrt_series!(res, res, length(res))
    return res
end

# Use that d/dt atan(y(t) / x(t)) = (y' * x - y * x') / (x^2 + y^2)
function Base.atan(y::ArbSeries, x::ArbSeries)
    res = zero(x)

    # res = y' * x - y * x'
    tmp1 = Arblib.derivative(y)
    Arblib.mullow!(tmp1, tmp1, x, length(tmp1))
    tmp2 = Arblib.derivative(x)
    Arblib.mullow!(tmp2, tmp2, y, length(tmp2))
    Arblib.sub_series!(res, tmp1, tmp2, length(res) - 1)

    # tmp1 = x^2 + y^2
    Arblib.mullow!(tmp1, x, x, length(tmp1))
    Arblib.mullow!(tmp2, y, y, length(tmp2))
    Arblib.add_series!(tmp1, tmp1, tmp2, length(tmp1))

    Arblib.div_series!(res, res, tmp1, length(res) - 1)
    Arblib.integral!(res, res)

    res[0] = atan(Arblib.ref(y, 0), Arblib.ref(x, 0))

    return res
end
