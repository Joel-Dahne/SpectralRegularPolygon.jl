function eigenvalue_lower_upper_estimate(u::Eigenfunction{T}, λ::T) where {T}
    m = maximum_boundary_estimate(u, λ)
    n = norm_estimate(u, λ)

    μ = sqrt(area(u.domain)) * m / n

    lower = λ / (1 + μ)
    upper = λ / (1 - μ)

    return lower, upper
end


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
