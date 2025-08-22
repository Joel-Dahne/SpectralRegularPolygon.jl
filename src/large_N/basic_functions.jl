exppii(z::Arblib.AcbOrRef) = Arblib.exp_pi_i!(zero(z), z)
exppii(z::Arblib.ArbOrRef) = exppii(Acb(z))

"""
    abspow!(res, x, y)

Inplace version of [`abspow`](@ref).
"""
function abspow!(res::Arb, x::Arblib.ArbOrRef, y::Arb)
    iszero(y) && return Arblib.one!(res)

    if iszero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(res)
        Arblib.ispositive(y) && return Arblib.zero!(res)
        return Arblib.unit_interval!(res)
    end

    if Arblib.contains_zero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(res)
        upper = abs_ubound(Arb, x) # One extra allocation
        Arblib.pow!(upper, upper, y)
        Arblib.zero!(res)
        return Arblib.union!(res, res, upper)
    end

    if res === y
        # In this case we need an extra allocation to not overwrite y
        y = copy(y)
    end
    Arblib.abs!(res, x)
    return Arblib.pow!(res, res, y)
end

"""
    abspow(x::Arb, y::Arb)

Compute `abs(x)^y` in a way that works if `x` overlaps with zero.
"""
abspow(x::Arb, y::Arb) = abspow!(zero(x), x, y)

"""
    pow(z::Acb, y::Arb)

Compute `z^y` in a way that works if `z` overlaps with zero.

For `z` overlapping zero it uses that the result is bounded in
absolute value by `abs(z)^y`, which can be computed using
[`abspow`](@ref).
"""
function pow(z::Acb, y::Arb)
    if Arblib.contains_zero(z)
        return Arblib.add_error!(zero(z), abspow(abs_ubound(Arb, z), y))
    else
        return z^y
    end
end

"""
    logabspow(x::Arb, m::Integer, y::Union{Arb,Integer})

Compute `log(abs(x))^m * abs(x)^y` in a way that works for `x`
overlapping zero.

For `x` overlapping zero it makes use of monotonicity in `x`. For `x >
0`, `y != 0` and `m != 0` the critical points of `log(x)^m * x^y` are
at `x = exp(-m / y)` and possibly `x = 1` (depending on the value of
`m`). At `x = exp(-m / y)` the value is given by
```
log(exp(-m / y))^m * exp(-m / y)^y = (-m / y)^m * exp(-m)
```
and at `x = 1` it is zero for `m > 0` and non-finite for `m < 0`.
"""
function logabspow(x::Arb, m::Integer, y::Arb)
    iszero(m) && return abspow(x, y)

    if Arblib.contains_zero(x)
        if Arblib.ispositive(y)
            iszero(x) && return zero(x)

            # Evaluate at endpoints of x
            xᵤ = abs_ubound(Arb, x)
            res = Arblib.union(zero(x), log(xᵤ)^m * xᵤ^y)

            # Check if critical points are contained in x, if so
            # evaluate on it
            # First critical point
            critical_point = exp(-m / y)
            if Arblib.overlaps(x, critical_point)
                res = Arblib.union(res, (-m / y)^m * exp(Arb(-m)))
            end
            # Second critical point
            if Arblib.contains(x, 1)
                # If m > 0 then the value is zero, which is already
                # included in res. Otherwise the value is non-finite.
                m < 0 && indeterminate!(res)
            end

            return res
        elseif iszero(y) && m < 0
            iszero(x) && return zero(x)

            xᵤ = abs_ubound(Arb, x)

            # log(1) = 0 so we get an indeterminate value there
            xᵤ < 1 || return indeterminate(x)

            # Monotone for 0 < x < 1, evaluate on endpoints
            return Arblib.union(zero(x), log(xᵤ)^i)
        else
            # Non-finite at x = 0
            return indeterminate(x)
        end
    end

    return log(abs(x))^m * abspow(x, y)
end

logabspow(x::Arb, m::Integer, y::Integer) = logabspow(x, m, Arb(y))

"""
    logpow(z::Acb, m::Integer, y::Union{Arb,Integer})

Compute `log(z)^m * z^y` in a way that works when `z` contains zero.

For `z` overlapping zero it uses that
```
log(z)^m = (log(abs(z)) + im * angle(z))^m
```
together with the binomial theorem. This gives a sum with terms of the
form
```
binomial(m, n) * log(abs(z))^n * (im * angle(z))^(m - n)
```
Putting `z^y` into the sum we then bound the factor
```
log(abs(z))^n * z^y
```
by writing it as
```
log(abs(z))^n * z^y = log(abs(z))^n * abs(z)^y * exp(im * angle(z))
```
and using [`logabspow`](@ref).
"""
function logpow(z::Acb, m::Integer, y::Arb)
    m >= 0 || throw(ArgumentError("only supports m >= 0"))

    m == 0 && return pow(z, y)

    if Arblib.contains_zero(z)
        return sum(0:m) do n
            binomial(m, n) *
            logabspow(abs(z), n, y) *
            exp(Acb(0, Base.angle(z))) *
            Acb(0, Base.angle(z))^(m - n)
        end
    else
        return log(z)^m * z^y
    end
end

logpow(z::Acb, i::Integer, y::Integer) = logpow(z, i, Arb(y))
