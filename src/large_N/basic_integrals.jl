"""
    integral_log(m::Int, a::Arb)

Compute the integral of `log(x)^m` from `0` to `a`.

It uses that the integral is given by
```
(-1)^m * a * sum(0:m) do n
    (-1)^n * factorial(m) / factorial(n) * log(a)^n
end
```
"""
function integral_log(m::Int, a::Arb)
    0 < a < 1 || throw(ArgumentError("only supports 0 < a < 1"))

    # Write sum as a polynomial for efficient evaluation
    p = ArbPoly([(-1)^n * factorial(m) // factorial(n) for n = 0:m])

    return (-1)^m * a * p(log(a))
end

"""
    integral_log_1mtz(z::Acb, m::Int, b::Arb)

Compute the integral of `log(1 - t * z)^m` from `b` to `1`.

It uses that the primitive function is given by
```
(-1)^m * factorial(m) * t - (-1)^m * (1 - t * z) / z * sum(1:m) do n
    (-1)^n * factorial(m) / factorial(n) * log(1 - t * z)^n
end
```
Which we can also write as
```
(-1)^m * factorial(m) * (
    t - 1 / z * sum(1:m) do n
        (-1)^n * log(1 - t * z)^n * (1 - t * z) / factorial(n)
    end
)
```
and use [`logpow`](@ref) to evaluate `log(1 - t * z)^n * (1 - t * z)`.
"""
function integral_log_1mtz(z::Arblib.AcbOrRef, m::Integer, b::Arb)
    0 < b < 1 || throw(ArgumentError("only supports 0 < b < 1"))
    m >= 0 || throw(ArgumentError("only supports m >= 0"))

    primitive(t) =
        (-1)^m *
        factorial(m) *
        (t - 1 / z * sum(1:m, init = zero(z)) do n
            (-1)^n * logpow(1 - t * z, n, 1) / factorial(n)
        end)

    return primitive(one(b)) - primitive(b)
end

"""
    integral_logpow_1mtz(z::Acb, m::Int, y::Arb, b::Arb)

Compute the integral of `log(1 - t * z)^m * (1 - t * z)^y` from `b` to
`1`.

It uses that the primitive function is given by
```
(-1)^(m + 1) / ((1 + y)^(m + 1) * z) * (1 - t * z)^(1 + y) * sum(0:m) do n
    (-1)^n * factorial(m) / factorial(n) * (1 + y)^n * log(1 - t * z)^n
end
```
Which we can also write as
```
(-1)^(m + 1) * factorial(m) / z * sum(0:m) do n
    (-1)^n  * (1 + y)^(n - m - 1) * log(1 - t * z)^n * (1 - t * z)^(1 + y) / factorial(n)
end
```
and use [`logpow`](@ref) to evaluate `log(1 - t * z)^n * (1 - t * z)^(1 + y)`.
"""
function integral_logpow_1mtz(z::Acb, m::Int, y::Arb, b::Arb)
    0 < b < 1 || throw(ArgumentError("only supports 0 < b < 1"))
    m >= 0 || throw(ArgumentError("only supports m >= 0"))

    # Primitive function
    primitive(t) =
        (-1)^(m + 1) * factorial(m) / z * sum(0:m) do n
            (-1)^n * (1 + y)^(n - m - 1) * logpow(1 - t * z, n, 1 + y) / factorial(n)
        end

    return primitive(one(b)) - primitive(b)
end
