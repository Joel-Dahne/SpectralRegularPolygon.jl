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

"""
    iswide(x; cutoff = 10)

Return true if `x` is wide, in the meaning that the effective relative
accuracy of `x` measured in bits is more than `cutoff` lower than it's
precision. For `x` not of type `Arb` or `Acb` this always return
`false`. For `x` of type `ArbSeries` or `AcbSeries` it checks the
first coefficient.
"""
iswide(x::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}; cutoff = 10) =
    Arblib.rel_accuracy_bits(x) < precision(x) - cutoff
iswide(x::Union{ArbSeries,AcbSeries}; cutoff = 10) = iswide(Arblib.ref(x, 0); cutoff)
iswide(::Number; cutoff = 10) = false

"""
    <<(p::Union{ArbSeries,AcbSeries}, n::Integer)

Return `p` divided by `x^n`, updating the degree accordingly

It throws an error if the lower order coefficients are not all exactly
equal to zero.

Note that the naming is different from Arb where division by `x^n` is
referred to as right shift whereas here we call it a left shift.
"""
function Base.:(<<)(p::T, n::Integer) where {T<:Union{ArbSeries,AcbSeries}}
    n >= 0 || throw(ArgumentError("n needs to be non-negative, got $n"))
    for i = 0:(n-1)
        iszero(Arblib.ref(p, i)) ||
            throw(ArgumentError("coefficient $i not equal to zero, got $(p[i])"))
    end
    return Arblib.shift_right!(T(degree = Arblib.degree(p) - n, prec = precision(p)), p, n)
end

"""
    >>(p::Union{ArbSeries,AcbSeries}, n::Integer)

Return `p` multiplied by `x^n`, updating the degree accordingly.

Note that the naming is different from Arb where multiplication by
`x^n` is referred to as left shift whereas here we call it a right
shift.
"""
function Base.:(>>)(p::T, n::Integer) where {T<:Union{ArbSeries,AcbSeries}}
    n >= 0 || throw(ArgumentError("n needs to be non-negative, got $n"))
    return Arblib.shift_left!(T(degree = Arblib.degree(p) + n, prec = precision(p)), p, n)
end
