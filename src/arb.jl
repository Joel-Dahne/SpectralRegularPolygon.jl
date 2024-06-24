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
