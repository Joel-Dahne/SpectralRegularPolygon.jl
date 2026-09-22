# Functions for handling V_l(z) from the paper.

"""
    V(l::Int, z::Arblib.AcbOrRef; analytic::Bool = false)

Compute an enclosure of `V_l(z)` from Equation REF(11) in the paper.

If `analytic` is true, then return an indeterminate value if `z`
overlaps a branch cut of the function.
"""
V(l::Int, z::Arblib.AcbOrRef; analytic::Bool = false) =
    if l == 1
        V_1(z; analytic)
    elseif l == 2
        V_2(z; analytic)
    elseif l == 3
        V_3(z; analytic)
    elseif l == 4
        V_4(z; analytic)
    end

"""
    V_log_bound_coefficients(l::Int)

Return coefficients `Cs` such that
```
abs(V(l, t)) <= sum(j -> C[j] * abs(log(1 - t))^j / factorial(j), 1:l)
```
for `abs(z) < 1`.

This is based on Lemma REF(A.2) in the paper.
"""
V_log_bound_coefficients(l::Int) =
    if l == 1
        V_1_log_bound_coefficients()
    elseif l == 2
        V_2_log_bound_coefficients()
    elseif l == 3
        V_3_log_bound_coefficients()
    elseif l == 4
        V_4_log_bound_coefficients()
    end

"""
    V_div_z_bound(l::Int, a::Arb)

Return a constant `D` such that
```
abs(V(l, z) / z) <= D
```
for all `|z| <= a < 1`.

This is based on Lemma REF(A.2) in the paper.
"""
function V_div_z_bound(l::Int, a::Arb)
    C = V_log_bound_coefficients(l)
    return sum(j -> C[j] * abs(log(1 - a))^j / (factorial(j) * a), 1:l)
end

function V_1(z; analytic::Bool = false)
    # Check analyticity when requested. It has a branch cut for z in
    # [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(z)) &&
       Arblib.contains_nonnegative(Arblib.realref(z) - 1)

        return indeterminate(z)
    end

    return 2polylog(1, z)
end

V_1_log_bound_coefficients() = Arb[2]

function V_2(z; analytic::Bool = false)
    # Check analyticity when requested. It has a branch cut for z in
    # [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(z)) &&
       Arblib.contains_nonnegative(Arblib.realref(z) - 1)

        return indeterminate(z)
    end

    λ = λ_disc()
    return (λ / 2 - 2) * polylog(2, z) + 4polylog_1_1(z)
end

function V_2_log_bound_coefficients()
    λ = λ_disc()
    return Arb[abs(λ / 2 - 2), 4]
end

function V_3(z; analytic::Bool = false)
    # Check analyticity when requested. It has a branch cut for z in
    # [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(z)) &&
       Arblib.contains_nonnegative(Arblib.realref(z) - 1)

        return indeterminate(z)
    end

    λ = λ_disc()
    return (λ^2 / 16 - λ + 2) * polylog(3, z) +
           (3λ - 12) * polylog_1_2(z) +
           (λ - 4) * polylog_2_1(z) +
           8polylog_1_1_1(z)
end

function V_3_log_bound_coefficients()
    λ = λ_disc()
    return Arb[abs(λ^2 / 16 - λ + 2), abs(3λ-12)+abs(λ-4), 8]
end

function V_4(z; analytic::Bool = false)
    # Check analyticity when requested. It has a branch cut for z in
    # [1, Inf].
    if analytic &&
       Arblib.contains_zero(Arblib.imagref(z)) &&
       Arblib.contains_nonnegative(Arblib.realref(z) - 1)

        return indeterminate(z)
    end

    λ = λ_disc()
    return (λ^3 / 192 - λ^2 / 8 - λ / 2 - 2) * polylog(4, z) +
           (λ^2 / 8 - 2λ + 4) * polylog_3_1(z) +
           (λ^2 / 4 - 4λ + 12) * polylog_2_2(z) +
           (5λ^2 / 8 - 8λ + 28) * polylog_1_3(z) +
           (2λ - 8) * polylog_2_1_1(z) +
           (6λ - 24) * polylog_1_2_1(z) +
           (14λ - 56) * polylog_1_1_2(z) +
           16polylog_1_1_1_1(z) +
           2λ * zeta(Arb(3)) * polylog(1, z)
end

function V_4_log_bound_coefficients()
    λ = λ_disc()
    return Arb[
        abs(λ^3/192-λ^2/8-λ/2-2)+abs(2λ*zeta(Arb(3))),
        abs(λ^2/8-2λ+4)+abs(λ^2/4-4λ+12)+abs(5λ^2/8-8λ+28),
        abs(2λ-8)+abs(6λ-24)+abs(14λ-56),
        16,
    ]
end
