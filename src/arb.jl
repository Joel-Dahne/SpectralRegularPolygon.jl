# FIXME: Don't overload a Base function
function Base.atan(y::ArbSeries, x::ArbSeries)
    res = atan(y / x)
    res[0] = atan(y[0], x[0])
    return res
end
