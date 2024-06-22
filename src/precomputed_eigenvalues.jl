function load_precomputed_eigenvalue(T, N::Integer)
    filename = "lowest_eigenvalues.txt"

    types = [Int, Int, Int, String]

    df = CSV.read(filename, DataFrames.DataFrame; types)

    row_index = searchsortedfirst(df[:, 1], N)

    if T == Arb
        return Arb(df.lambda[row_index])
    else
        return parse(T, df.lambda[row_index])
    end
end
