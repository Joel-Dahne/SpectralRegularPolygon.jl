abstract type AbstractVertexExpansion{T} end
abstract type AbstractInteriorExpansion{T} end

struct VertexExpansion{T} <: AbstractVertexExpansion{T}
    θ::T
    coefficients::Vector{T}
end

struct InteriorExpansion{T} <: AbstractInteriorExpansion{T}
    coefficients::Vector{T}
    symmetry::Int
end

struct Eigenfunction{T,U<:AbstractVertexExpansion{T},V<:AbstractInteriorExpansion{T}}
    domain::RegularPolygon{T}
    vertex_expansion::U
    interior_expansion::V
end

function Eigenfunction(
    domain::RegularPolygon{T},
    vertex_expansion::U,
    interior_expansion::V,
) where {T,U,V}
    u = Eigenfunction{T,U,V}(domain, vertex_expansion, interior_expansion)

    return u
end

function Eigenfunction(domain::RegularPolygon{T}; lightning::Bool = false) where {T}
    vertex_expansion = if lightning
        LightningExpansion{T}(T(1.0), T(2.5), angle(domain, 1), true, T[])
    else
        VertexExpansion{T}(angle(domain, 1), T[])
    end

    u = Eigenfunction(domain, vertex_expansion, InteriorExpansion{T}(T[], domain.N))

    return u
end

num_coefficients(u::Eigenfunction) =
    length(u.vertex_expansion.coefficients) + length(u.interior_expansion.coefficients)

function Base.show(io::IO, ::MIME"text/plain", u::Eigenfunction{T}) where {T}
    n = num_coefficients(u)
    print(io, "$(typeof(u)) with $n coefficients")
end

function set_coefficients!(u::Eigenfunction{T}, coefficients::Vector{T}) where {T}
    # Distribute the coefficients alternating between the vertex
    # expansion and the interior expansion.
    # IMPROVE: Allow other splitting
    copy!(u.vertex_expansion.coefficients, coefficients[1:2:end])
    copy!(u.interior_expansion.coefficients, coefficients[2:2:end])

    return u
end

"""
    normalise_sign!(u::Eigenfunction{T})

Normalise the sign of the eigenfunction so that it is positive at the
center of the domain.
"""
function normalise_sign!(u::Eigenfunction{T}, λ::T) where {T}
    if u(center(u.domain), λ) < 0
        u.vertex_expansion.coefficients .*= -1
        u.interior_expansion.coefficients .*= -1
    end

    return u
end

function (v::VertexExpansion{T})(p::Polar{T}, λ::T, ks::UnitRange{Int}) where {T}
    r_sqrt_λ = p.r * sqrt(λ)
    π_div_θ = π / v.θ

    return map(ks) do k
        ν = (1 + 2(k - 1)) * π_div_θ
        besselj(ν, r_sqrt_λ) * sin(ν * p.φ)
    end
end

function (v::VertexExpansion{T})(p::Polar, λ::T) where {T}
    r_sqrt_λ = p.r * sqrt(λ)
    π_div_θ = π / v.θ # IMPROVE: Issues if this is exact integer

    return sum(eachindex(v.coefficients), init = zero(r_sqrt_λ)) do k
        ν = (1 + 2(k - 1)) * π_div_θ
        v.coefficients[k] * besselj(ν, r_sqrt_λ) * sin(ν * p.φ)
    end
end

function (v::InteriorExpansion{T})(p::Polar{T}, λ::T, ks::UnitRange{Int}) where {T}
    r_sqrt_λ = p.r * sqrt(λ)

    return map(ks) do k
        if k == 1
            ν = zero(T)
            term = besselj(ν, r_sqrt_λ)
        else
            ν = convert(T, v.symmetry * (k - 1))
            term = besselj(ν, r_sqrt_λ) * cos(ν * p.φ)
        end

        term
    end
end

function (v::InteriorExpansion{T})(p::Polar, λ::T) where {T}
    r_sqrt_λ = p.r * sqrt(λ)

    return sum(eachindex(v.coefficients), init = zero(r_sqrt_λ)) do k
        if k == 1
            ν = zero(T)
            term = besselj(ν, r_sqrt_λ)
        else
            ν = convert(T, v.symmetry * (k - 1))
            term = besselj(ν, r_sqrt_λ) * cos(ν * p.φ)
        end

        v.coefficients[k] * term
    end
end

# Return which vertices for which evaluation of the given boundary
# point is non-zero
active_vertices(
    v::VertexExpansion{T},
    domain::RegularPolygon{T},
    xy::BoundaryPoint2,
) where {T} = mod1.(xy.boundary .+ (2:domain.N-1), domain.N)

active_vertices(
    v::LightningExpansion{T},
    domain::RegularPolygon{T},
    xy::BoundaryPoint2,
) where {T} = 1:domain.N

function (u::Eigenfunction{T})(
    xy::Union{Point2{T},BoundaryPoint2{T}},
    λ::T,
    n::Integer,
) where {T}
    # Distribute the coefficients alternating between the vertex
    # expansion and the interior expansion.
    n_per = n ÷ 2
    ks_1 = 1:(n_per+(n%2==1))
    ks_2 = 1:n_per

    if xy isa Point2
        res_vertices = sum(1:u.domain.N) do i
            u.vertex_expansion(polar_vertex(u.domain, xy, i), λ, ks_1)
        end

        res_interior = u.interior_expansion(polar_center(u.domain, xy), λ, ks_2)
    else
        res_vertices = sum(active_vertices(u.vertex_expansion, u.domain, xy)) do i
            u.vertex_expansion(polar_vertex(u.domain, xy.position, i), λ, ks_1)
        end

        res_interior = u.interior_expansion(polar_center(u.domain, xy.position), λ, ks_2)
    end

    ress = (res_vertices, res_interior)

    # Collect the computed values in an inverse round robin fashion.
    res = Vector{T}(undef, n)
    for k = 1:n
        res[k] = ress[mod1(k, 2)][1+(k-1)÷2]
    end

    return res
end

function (u::Eigenfunction{T})(xy::Point2, λ::T) where {T}
    res_vertices = sum(1:u.domain.N) do i
        u.vertex_expansion(polar_vertex(u.domain, xy, i), λ)
    end

    res_interior = u.interior_expansion(polar_center(u.domain, xy), λ)

    return res_vertices + res_interior
end

function (u::Eigenfunction{T})(xy::BoundaryPoint2, λ::T) where {T}
    res_vertices = sum(active_vertices(u.vertex_expansion, u.domain, xy)) do i
        u.vertex_expansion(polar_vertex(u.domain, xy.position, i), λ)
    end

    res_interior = u.interior_expansion(polar_center(u.domain, xy.position), λ)

    return res_vertices + res_interior
end

(u::Eigenfunction{T})(_::Missing, _::T) where {T} = missing
