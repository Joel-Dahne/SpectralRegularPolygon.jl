"""
    sigma_matrix(u, λ, N; num_boundary, num_interior)

Compute the matrix which is used for computing ``σ(λ)`` in the MPS.

The matrix is of size `(num_boundary + num_interior, N)`.

The first `num_boundary` rows consists of the first `N` basis
functions of `u` evaluated on boundary points of the domain. The
remaining `num_interior` rows corresponds to the basis functions
evaluated on interior points.
"""
function sigma_matrix(
    u::Eigenfunction{T},
    λ::T,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
) where {T}
    A = Matrix{T}(undef, num_boundary + num_interior, N)

    boundary = boundary_points_symmetry(u.domain, num_boundary)
    for i = 1:num_boundary
        A[i, :] = u(boundary[i], λ, N)
    end

    interior = interior_points_random(u.domain, num_interior)
    for i = 1:num_interior
        A[num_boundary+i, :] = u(interior[i], λ, N)
    end

    return A
end

"""
    sigma(u, λ, N)

Compute ``σ(λ)``.

# Numerical stability
In some cases when using the lightning expansions there are issues
with the numerical stability in the computations. For example for
```
u = Eigenfunction(RegularPolygon{T}(12), lightning = true)
```
with `N = 8`.

It seems like these numerical instabilities come from the
QR-decomposition.

If we compute the QR-decomposition at a fixed precision (say 2048
bits) and then compute the SVD-values at lower precisions (say
`64:64:2048` (we first truncate to the lower precision)), we get that
the SVD-values are stable.

If we instead compute the QR-decomposition at different precisions
(say `64:64:2048`), and then compute the SVD-values at high precision
(say 2048) we get that they vary a lot.

"""
function sigma(
    u::Eigenfunction{T},
    λ::T,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
    qr_eltype = ifelse(T == Arb, BigFloat, T),
) where {T}
    # Compute the matrix A
    A = convert(Matrix{qr_eltype}, sigma_matrix(u, λ, N; num_boundary, num_interior))

    # Compute a QR factorization of A
    Q = Matrix(LinearAlgebra.qr(A).Q)

    # Compute the smallest singular value of the top part of Q,
    # corresponding to the boundary points
    σ = LinearAlgebra.svdvals(Q[1:num_boundary, :])[end]

    return convert(T, σ)
end

"""
    sigma!(u, λ, N)

Compute ``σ(λ)`` and set the coefficients of `u` to the corresponding
ones.
"""
function sigma!(
    u::Eigenfunction{T},
    λ::T,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
    qr_eltype = ifelse(T == Arb, BigFloat, T),
    normalise::Bool = true,
) where {T}
    # Compute the matrix A
    A = convert(Matrix{qr_eltype}, sigma_matrix(u, λ, N; num_boundary, num_interior))

    # Compute a QR factorization of A
    q = LinearAlgebra.qr(A)
    Q = Matrix(q.Q)

    # Compute the svd of the top part of Q, corresponding to the
    # boundary points
    s = LinearAlgebra.svd(Q[1:num_boundary, :])

    # Extract smallest singular value and corresponding vector
    σ = s.S[end]
    v = s.V[:, end]

    # Compute the coefficients
    coefficients = try
        q \ (Q * v)
    catch
        @warn "Failed computing q \\ (Q * v)"
        zeros(qr_eltype, N)
    end

    set_coefficients!(u, convert(Vector{T}, coefficients))

    if normalise
        normalise_leading_interior!(u)
    end

    return convert(T, σ)
end
