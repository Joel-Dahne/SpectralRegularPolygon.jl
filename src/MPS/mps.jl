"""
    mps(
        u::Eigenfunction{T},
        a::T,
        b::T,
        N::Integer;
        num_boundary::Integer = 2N,
        num_interior::Integer = 2N,
        qr_eltype = ifelse(T == Arb, BigFloat, T),
        verbose = false,
)

Find `λ` between `a` and `b` that minimizes `sigma(u, λ, N)`.
"""
function mps(
    u::Eigenfunction{T},
    a::T,
    b::T,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
    qr_eltype = ifelse(T == Arb, BigFloat, T),
    verbose = false,
) where {T}
    res = Optim.optimize(a, b, show_trace = verbose) do λ
        sigma(u, λ, N; num_boundary, num_interior, qr_eltype)
    end

    Optim.converged(res) || @warn "Optimization failed"

    return Optim.minimizer(res)
end

"""
    mps(
        u::Eigenfunction{T},
        λ₀::T,
        N::Integer;
        num_boundary::Integer = 2N,
        num_interior::Integer = 2N,
        qr_eltype = ifelse(T == Arb, BigFloat, T),
        verbose = false,
)

Find `λ` minimizes `sigma(u, λ, N)`, starting from the initial guess `λ₀`.

Note that this is significantly slower than the version that takes `a`
and `b`. It can however be convinient if one has a decent
approximation.
"""
function mps(
    u::Eigenfunction{T},
    λ₀::T,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
    qr_eltype = ifelse(T == Arb, BigFloat, T),
    verbose = false,
) where {T}
    # TODO: This might not be a good idea
    res = Optim.optimize([λ₀], LBFGS()) do λ
        sigma(u, only(λ), N; num_boundary, num_interior, qr_eltype)
    end

    Optim.converged(res) || @warn "Optimization failed"

    return only(Optim.minimizer(res))
end

"""
    mps!(
        u::Eigenfunction{T},
        a::T,
        b::T,
        N::Integer;
        num_boundary::Integer = 2N,
        num_interior::Integer = 2N,
        qr_eltype = ifelse(T == Arb, BigFloat, T),
        verbose = false,
    )

Like [`mps`](@ref) but also sets the coefficients of `u` to the ones
corresponding to `λ`.
"""
function mps!(
    u::Eigenfunction{T},
    a::T,
    b::T,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
    qr_eltype = ifelse(T == Arb, BigFloat, T),
    verbose = false,
) where {T}
    λ = mps(u, a, b, N; num_boundary, num_interior, qr_eltype, verbose)

    sigma!(u, λ, N; num_boundary, num_interior, qr_eltype)

    return λ
end

function mps!(
    u::Eigenfunction{T},
    λ₀::T,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
    qr_eltype = ifelse(T == Arb, BigFloat, T),
    verbose = false,
) where {T}
    λ = mps(u, λ₀, N; num_boundary, num_interior, qr_eltype, verbose)

    sigma!(u, λ, N; num_boundary, num_interior, qr_eltype)

    return λ
end
