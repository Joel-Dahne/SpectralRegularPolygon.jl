"""
    norm_estimate(u::Eigenfunction{T}, λ::T, n::Integer)

Compute an estimate of the norm of `u` by evaluating it on `n` random
points.
"""
function norm_estimate(u::Eigenfunction{T}, λ::T, n::Integer = 100) where {T}
    interior = interior_points_random(u.domain, n)
    sqrt(area(u.domain) * sum(abs(u(xy, λ))^2 for xy in interior) / length(interior))
end

"""
    norm_lower_estimate(u::Eigenfunction{T}, λ::T, scaling::T = T(0.65), n::Integer = 100)

Compute an estimate of a lower bound of the norm of `u`.

The lower bound is given by considering the norm of `u` in a scaled
down version of the original domain. A lower bound of the norm in that
subdomain is computed by lower bounding `u` on its boundary and using
the maximum principle.

The lower bound of `u` on the boundary of the subdomain is estimated
by evaluating it on `n` values.

The default `scaling = 0.65` has been chosen to be reasonably close to
the optimal value for the circle.
"""
function norm_lower_estimate(
    u::Eigenfunction{T},
    λ::T,
    scaling::T = T(0.65),
    n::Integer = 100,
) where {T}
    u_lower_bound = maximum(boundary_points_symmetry(u.domain, n), init = zero(T)) do xy
        abs(u(scaling * xy.position, λ))
    end

    return area(u.domain) * scaling^2 * u_lower_bound
end

"""
    norm_lower_enclosure(u::Eigenfunction{T}, λ::T, scaling::T = T(0.65), n::Integer = 100)

Compute an enclosure of a lower bound of the norm of `u`.

The lower bound is given by considering the norm of `u` in a scaled
down version of the original domain. A lower bound of the norm in that
subdomain is computed by lower bounding `u` on its boundary and using
the maximum principle.

The default `scaling = 0.65` has been chosen to be reasonably close to
the optimal value for the circle.
"""
function norm_lower_enclosure(
    u::Eigenfunction{T},
    λ::T,
    scaling::T = T(0.65);
    degree::Integer = num_coefficients(u),
    rtol = 1e-3,
    maxevals::Integer = 100,
    depth::Integer = 10,
    threaded::Bool = true,
    verbose::Bool = false,
) where {T}
    u_lower_bound = ArbExtras.maximum_enclosure(
        t -> u(scaling * boundary_parameterized_symmetry(u.domain, t).position, λ),
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

    return area(u.domain) * scaling^2 * u_lower_bound
end
