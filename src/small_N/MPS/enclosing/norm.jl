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
subdomain is computed by lower bounding `abs(u)` on its boundary and
using the maximum principle.

For the maximum principle to be applicable we need `u` to be non-zero
on the boundary as well as in the subdomain. If the lower bound of
`abs(u)` on the boundary is non-zero then `u` is non-zero on the
boundary. To ensure that it is non-zero in the interior of the
subdomain we make use of the Faber-Krahn inequality, however this part
is not checked in this method since it is only intended to compute an
estimate.

The default `scaling = 0.65` has been chosen to be reasonably close to
the optimal value for the circle.
"""
function norm_lower_estimate(
    u::Eigenfunction{T},
    λ::T,
    scaling::T = T(0.65),
    n::Integer = 100,
) where {T}
    u_lower_bound = minimum(boundary_points_symmetry(u.domain, n), init = zero(T)) do xy
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

For the maximum principle to be applicable we need `u` to be non-zero
on the boundary as well as in the subdomain. If the lower bound of
`abs(u)` on the boundary is non-zero then `u` is non-zero on the
boundary. To ensure that it is non-zero in the interior of the
subdomain we make use of the Faber-Krahn inequality. We are guaranteed
that `u` is non-zero in the domain if `λ` is smaller than the first
eigenvalue of the disc with the same area as the subdomain.

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
    # Check Faber-Krahn condition
    λ_unit_disc = let
        # The first zero of besselj0 is approximately 2.4048, and is
        # unique on the interval [0, 2.5]. To compute an enclosure we
        # verify that there is a unique root on the interval [0, 2.5]
        # and enclose it.
        zeros, flags = ArbExtras.isolate_roots(besselj0, Arf(0), Arf(2.5))
        # Verify uniqueness
        only(flags) || error("could not prove unique root on [0, 2.5] for besselj0")

        ArbExtras.refine_root(besselj0, Arb(only(zeros)))^2
    end
    λ_scaled_disc = λ_unit_disc / scaling^2

    if !(λ < λ_scaled_disc)
        error("Faber-Krahn condition not satisfied")
    end

    # Lower bound on boundary
    u_lower_bound = ArbExtras.minimum_enclosure(
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

    # Note that if u_lower_bound contains zero then the returned value
    # also contains zero, so this is a lower bound even if all
    # conditions are not met.
    return area(u.domain) * scaling^2 * u_lower_bound
end
