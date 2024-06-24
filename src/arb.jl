# FIXME: Don't overload a Base function
function Base.atan(y::ArbSeries, x::ArbSeries)
    # Use that d/dt atan(y(t) / x(t)) =
    # (y' * x - y * x') / (x^2 + y^2)
    res =
        Arblib.integral((Arblib.derivative(y) * x - y * Arblib.derivative(x)) / (x^2 + y^2))
    res[0] = atan(y[0], x[0])

    return res
end
