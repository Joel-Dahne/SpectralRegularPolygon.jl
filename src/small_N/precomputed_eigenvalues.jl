function get_eigenvalue_approximation(T, N::Integer)
    if N == 3
        return 4T(π) / sqrt(T(3))
    elseif N == 4
        return 2T(π)
    elseif 5 <= N <= 11
        return T(NaN)
    elseif 12 <= N <= 660
        return load_precomputed_eigenvalue(T, N)
    else
        return T(NaN)
    end
end

function load_precomputed_eigenvalue(T, N::Integer)
    filename = joinpath(dirname(pathof(@__MODULE__)), "../", "lowest_eigenvalues.txt")
    types = [Int, Int, Int, String]

    df = CSV.read(filename, DataFrames.DataFrame; types)

    row_index = searchsortedfirst(df[:, 1], N)

    if T == Arb
        return Arb(df.lambda[row_index])
    else
        return parse(T, df.lambda[row_index])
    end
end
