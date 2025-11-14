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
