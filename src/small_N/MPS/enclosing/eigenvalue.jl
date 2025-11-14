function eigenvalue_lower_upper_estimate(u::Eigenfunction{T}, λ::T) where {T}
    m = maximum_boundary_estimate(u, λ)
    n = norm_estimate(u, λ)

    μ = sqrt(area(u.domain)) * m / n

    lower = λ / (1 + μ)
    upper = λ / (1 - μ)

    return lower, upper
end

"""
    eigenvalue_enclosure(u::Eigenfunction{Arb}, λ::Arb; threaded = true, verbose = false)

Given the approximate eigenvalue `λ` and with associated approximate
eigenfunction `u`, compute an enclosure of true eigenvalue in the
neighborhood of `λ`.

The enclosure is given by computing `μ = sqrt(area(u.domain)) * m /
n`, where `m` is an upper bound of `u` in the boundary computed using
[`maximum_boundary_enclosure`](@ref) and `n` is a lower bound of the
norm of `u` computed using [`norm_lower_enclosure`](@ref). By Lemma
2.2 in the paper there is then an eigenvalue in the interval ``[λ / (1
+ μ), λ / (1 - μ)]``.
"""
function eigenvalue_enclosure(
    u::Eigenfunction{Arb},
    λ::Arb;
    threaded::Bool = true,
    verbose::Bool = false,
)
    verbose && @info "Computing maximum on boundary"
    m = maximum_boundary_enclosure(u, λ; threaded, verbose)
    verbose && @info "Computing lower bound of norm"
    n = norm_lower_enclosure(u, λ; threaded, verbose)

    μ = sqrt(area(u.domain)) * m / n

    lower = λ / (1 + μ)
    upper = λ / (1 - μ)

    return Arb((lower, upper))
end
