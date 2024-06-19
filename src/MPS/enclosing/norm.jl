"""
    norm_estimate(u::Eigenfunction{T}, λ::T, n::Integer)

Compute an estimate of the norm of `u` by evaluating it on `n` random
points.
"""
function norm_estimate(u::Eigenfunction{T}, λ::T, n::Integer = 100) where {T}
    interior = interior_points_random(u.domain, n)
    sqrt(area(u.domain) * sum(abs(u(xy, λ))^2 for xy in interior) / length(interior))
end
