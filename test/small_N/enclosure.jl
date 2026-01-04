@testset "enclosure" begin
    for N = 13:13:26
        domain = SRP.RegularPolygon{Arb}(N)
        u = SRP.Eigenfunction(domain)

        λ_approx = SRP.get_eigenvalue_approximation(Arb, N)

        perturbations = [-reverse(logrange(1e-10, 1e-4, 3)); 0; logrange(1e-10, 1e-4, 3)]
        Ms = 3:2:9

        λs = map(Iterators.product(perturbations, Ms)) do (perturbation, M)
            SRP.sigma!(u, λ_approx + perturbation, M)
            λ = SRP.eigenvalue_enclosure(u, λ_approx + perturbation)
        end

        # Compute intersection of all enclosures. If they do not all
        # intersect then this will throw an error.
        @test isfinite(foldl(Arblib.intersection, λs[:]))
    end
end
