"""
    maximum_boundary_estimate(u::Eigenfunction{T}, λ::T; num_points = 9num_coefficients(u))

Compute an estimate of the maximum of `abs(u)` on the boundary of the
domain.

The estimate is computed by evaluating it on `num_points` boundary
points.
"""
maximum_boundary_estimate(
    u::Eigenfunction{T},
    λ::T;
    num_points = 9num_coefficients(u),
) where {T} =
# FIXME: Take into account symmetries when taking boundary points
    maximum(boundary_points(u.domain, num_points), init = zero(T)) do xy
        abs(u(xy, λ))
    end
