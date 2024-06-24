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
    num_points = 8num_coefficients(u),
) where {T} =
    maximum(boundary_points_symmetry(u.domain, num_points), init = zero(T)) do xy
        abs(u(xy, λ))
    end

"""
    maximum_boundary_enclosure(u::Eigenfunction{Arb}, λ::Arb)

Compute an enclosure of the maximum of `abs(u)` on the boundary of the
domain.
"""
maximum_boundary_enclosure(
    u::Eigenfunction{Arb},
    λ::Arb;
    degree::Integer = 2num_coefficients(u),
    rtol = 1e-3,
    maxevals::Integer = 100,
    depth::Integer = 10,
    threaded::Bool = true,
    verbose::Bool = false,
) = ArbExtras.maximum_enclosure(
    t -> u(boundary_parameterized_symmetry(u.domain, t), λ),
    Arf(0),
    Arf(1),
    abs_value = true;
    degree,
    rtol,
    maxevals,
    depth,
    threaded,
    verbose,
)
