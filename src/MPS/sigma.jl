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
"""
function sigma(
    u::Eigenfunction{T},
    λ::T,
    N::Integer;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
) where {T}
    # Compute the matrix A
    A = sigma_matrix(u, λ, N; num_boundary, num_interior)

    # Compute a QR factorization of A
    Q = Matrix(LinearAlgebra.qr(A).Q)

    # Compute the smallest singular value of the top part of Q,
    # corresponding to the boundary points
    return LinearAlgebra.svdvals(Q[1:num_boundary, :])[end]
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
) where {T}
    # Compute the matrix A
    A = sigma_matrix(u, λ, N; num_boundary, num_interior)

    if T == Arb
        A = BigFloat.(A)
    end

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
        convert.(T, q \ (Q * v))
    catch
        @warn "Failed computing q \\ (Q * v)"
        zeros(T, N)
    end

    set_coefficients!(u, coefficients)

    return convert(T, σ)
end
