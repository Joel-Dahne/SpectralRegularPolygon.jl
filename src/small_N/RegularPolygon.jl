"""
    RegularPolygon{T<:Real}(N)

Represents a regular polygon of area `π` with `N` vertices.

Computations with the polygon return values of type `T`.
"""
struct RegularPolygon{T<:Real}
    N::Int

    function RegularPolygon{T}(N) where {T}
        N >= 3 || throw(DomainError(N, "N must be at least 3"))

        return new{T}(N)
    end
end

RegularPolygon(N) = RegularPolygon{Float64}(N)

angle(domain::RegularPolygon{T}, i::Integer) where {T} = π * T((domain.N - 2) // domain.N)

function vertex(domain::RegularPolygon{T}, i::Integer) where {T}
    s, c = sincospi(T(1 // domain.N))
    r = sqrt(π / (domain.N * s * c))

    s, c = sincospi(T(2 * (2(i - 1) + 1) // 2domain.N))

    return r * Point(c, s)
end

vertices(domain::RegularPolygon) = [vertex(domain, i) for i = 1:domain.N]

center(domain::RegularPolygon{T}) where {T} = Point(zero(T), zero(T))

# Should always be π with the way we currently choose the vertices
area(domain::RegularPolygon{T}) where {T} =
    domain.N // 4 * (2vertex(domain, 1)[2])^2 * cot(π / T(domain.N))

"""
    BoundaryPoint2{T}(position; boundary)

Represents a point `p` on a specified boundary of a [`RegularPolygon`](@ref).

The location of the point is at position `p.position` and it lies
along boundary `p.boundary` of the polygon.
"""
const BoundaryPoint2{T} = GeometryBasics.PointMeta{2,T,Point2{T},(:boundary,),Tuple{Int64}}

"""
    boundary_points(domain::RegularPolygon, i::Integer, n::Integer)

Return `n` points from boundary number `i`.
"""
function boundary_points(domain::RegularPolygon{T}, i::Integer, n::Integer) where {T}
    v = vertex(domain, i)
    w = vertex(domain, i + 1)

    points = Vector{BoundaryPoint2{T}}(undef, n)
    for j = 1:n
        # We currently only make use of Chebyshev spaced points, the
        # root exponential version is here for reference.
        if true
            # Take Chebyshev spaced points
            t = 1 - (cospi(T((2j - 1) // 2n)) + 1) / 2
        else
            # Take root exponentially spaced points
            m = (n + 1) // 2
            d = exp(-4(sqrt(T(m)) - sqrt(T(min(j, n - j + 1)))))
            t = if j < m
                d / 2
            elseif j > m
                1 - d / 2
            else
                T(0.5)
            end
        end

        points[j] = BoundaryPoint2{T}(v + t * (w - v), boundary = i)
    end

    return points
end

"""
    boundary_points(domain::RegularPolygon, n::Integer)

Return `n` points taken from all boundaries of the domain
"""
function boundary_points(domain::RegularPolygon, n::Integer)
    N = domain.N

    return reduce(vcat, [boundary_points(domain, i, n ÷ N + (n % N >= i)) for i = 1:N])
end

"""
    boundary_points_symmetry(domain::RegularPolygon, n::Integer)

Return `n` points from the boundary, taken in a way that takes into
account the symmetries.

For a regular polygon this means we take points along the upper part
of the boundary between the first and last vertices.
"""
boundary_points_symmetry(domain::RegularPolygon, n::Integer) =
    boundary_points(domain, domain.N, 2n - 1)[n:(2n-1)]

"""
    boundary_parameterized_symmetry(domain::RegularPolygon{T}, t::T)

Given a parameterization on the interval ``[0, 1]`` of the boundary of
the domain that takes into account the symmetry, return the point at
value `t`.

The parameterization is for the upper part of the boundary between
the first and last vertices.
"""
function boundary_parameterized_symmetry(domain::RegularPolygon, t)
    v = vertex(domain, 1)
    # Multiplication by one(t) is needed to support ArbSeries
    # correctly.
    position = Point(one(t) * v[1], t * v[2])
    return BoundaryPoint2{eltype(position)}(position, boundary = domain.N)
end

"""
    interior_points_random(domain::RegularPolygon, n::Integer; rng = Random.MersenneTwister(42))

Return `n` points randomly sampled from the interior of the domain
"""
function interior_points_random(
    domain::RegularPolygon{T},
    n::Integer;
    rng = Random.MersenneTwister(42),
) where {T}
    v = vertex(domain, 1)
    w = vertex(domain, domain.N)

    points = Vector{typeof(v)}(undef, n)
    for i = 1:n
        # Generate random point on the triangle formed between the
        # x-axis and the line between the origin and the first vertex.
        u1, u2 = rand(rng, Float64), rand(rng, Float64)
        if u1 + u2 > 1
            u1, u2 = 1 - u1, 1 - u2
        end

        triangle_point = u1 * v + u2 * w

        # Randomize rotation
        rotation = rand(rng, 0:(domain.N-1))

        s, c = sincospi(T(2 * rotation // domain.N))

        points[i] = Point2{T}(
            -s * triangle_point[2] + c * triangle_point[1],
            s * triangle_point[1] + c * triangle_point[2],
        )

    end

    return points
end

"""
    interior_points_grid(domain::RegularPolygon, n::Integer)

Return `n^2` points in a uniform square grid around the domain
together with boolean vector which indicates which points are inside
the domain.

The argument `dy` makes points count as inside if they are at most
`dy` distance away from the domain in the y-direction. This is
mostly used to play nicely with `heatmap`.
"""
function interior_points_grid(domain::RegularPolygon{T}, n::Integer; dr = 2 / n) where {T}
    n > 1 || throw(ArgumentError("requires n > 1, got n = $n"))

    x_min = vertex(domain, (domain.N + 1) ÷ 2)[1]
    x_max = vertex(domain, 1)[1]

    y_min = vertex(domain, (domain.N + 3) ÷ 4)[2]
    y_max = vertex(domain, (3domain.N + 3) ÷ 4)[2]

    max_norm = LinearAlgebra.norm(vertex(domain, 1)) + dr

    points = Vector{Union{Point2{T},Missing}}(undef, n^2)
    inside = similar(points, Bool)

    idx = 1
    for i = 1:n
        for j = 1:n
            points[idx] = Point(
                x_min + (i - 1) // (n - 1) * (x_max - x_min),
                y_min + (j - 1) // (n - 1) * (y_max - y_min),
            )
            # IMPROVE: Properly check if inside or not. Since this is
            # only used for plotting it is not essential to fix.
            inside[idx] = LinearAlgebra.norm(points[idx]) < max_norm

            idx += 1
        end
    end

    return points, inside
end

"""
    polar_vertex(domain::RegularPolygon, xy::Point2, i::Integer)

Convert from Cartesian coordinates to polar coordinates around vertex
`i` of the domain, with the angle taken to be zero along the edge
between vertex `i` and `i + 1`.
"""
function polar_vertex(domain::RegularPolygon{T}, xy::Point2, i::Integer) where {T}
    # Compute required rotation as a rational multiple of π
    angle = (domain.N - 2) // domain.N
    outer_angle = 2 - angle
    rotation = mod(-outer_angle + 1 // 2 - (i - 1) * (1 - angle), 2)

    return Polar(xy - vertex(domain, i), π * T(rotation))
end

"""
    polar_center(domain::RegularPolygon, xy::Point2)

Convert from Cartesian coordinates to polar coordinates around
`center(domain)`, with the angle taken to be zero in the direction of
the positive x-axis.
"""
polar_center(domain::RegularPolygon, xy::Point2) = Polar(xy - center(domain))
