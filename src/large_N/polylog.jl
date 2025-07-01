polylog(s::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}) = Arblib.polylog!(zero(z), s, z)

function polylog(s::Int, z::Union{ArbSeries,AcbSeries})
    z0 = z[0]
    res = zero(z)

    res[0] = polylog(s, z0)

    if length(z) > 1
        res[1] = polylog(s - 1, z0) / z0

        if length(z) > 2
            error("not implemented")
        end
    end

    return ArbExtras.compose_zero!(res, res, z)
end

"""
    polylog_disc(s::Int, z::Arblib.AcbOrRef)

Compute `polylog(s, z)` assuming that `abs(z) <= 1`. It is intended to
be used for `z` overlapping `1`, otherwise there are much more
efficient methods.

It uses the power series expansion at `z = 0` and bounds the tail by
```
sum(k -> 1 / k^s, N:Inf) = polygamma(1, N)
```

IMPROVE: The bound is very slowly converging. It might be sufficient
for what we need though.
"""
function polylog_disc(s::Int, z::Arblib.AcbOrRef)
    N = 10000

    res = sum(1:(N-1)) do k
        z^k / Arb(k)^s
    end

    tail = abs(Arblib.polygamma!(zero(z), one(z), Acb(N)))

    return Arblib.add_error!(res, tail)
end

function polylog_unsafe(s::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef})
    if Arblib.contains(z, one(z))
        # TODO: We here assume that this is only ever called with
        # z in the unit disc.
        return polylog_disc(s, z)
    else
        return polylog(s, z)
    end
end

polylog_unsafe(s::Int, z) = polylog(s, z)

lerch_phi(z::Acb, s::Int, a::Int) = Arblib.dirichlet_lerch_phi!(zero(z), z, Acb(s), Acb(a))

function S(n::Int, p::Int, z::Union{Arblib.ArbOrRef,Arblib.AcbOrRef,ArbSeries,AcbSeries})
    s = n + 1
    if p == 1
        return polylog_unsafe(s, z)
    elseif p == 2
        # FIXME: This is not rigorous and converges slowly for abs(z)
        # = 1. Should prefer to rewrite in terms of polylog using
        # recurrence.
        sum(1:1000) do k₁
            z^(k₁ + 1) / k₁ * lerch_phi(z, s, k₁ + 1)
        end
    elseif p == 3
        sum(1:100) do k₁
            sum((k₁+1):101) do k₂
                z^(k₂ + 1) / k₂ * lerch_phi(z, s, k₂ + 1)
            end / k₁
        end
    elseif p == 4
        sum(1:20) do k₁
            sum((k₁+1):21) do k₂
                sum((k₂+1):22) do k₃
                    z^(k₃ + 1) / k₃ * lerch_phi(z, s, k₃ + 1)
                end / k₂
            end / k₁
        end
    end
end

function S(n::Int, z)
    sum(1:(n-1)) do j
        (-1)^(j - 1) * 2^(n - j) * S(j, n - j, z)
    end
end
