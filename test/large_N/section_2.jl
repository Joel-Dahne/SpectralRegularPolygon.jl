@testset "section_2" begin
    @testset "F_N_model" begin
        N₀ = 64

        # Compute F_N using both F_N_model and its definition through
        # hypgeom2f1 and ensure that they agree.
        for z in Acb[1, SRP.exppii(Arb(1 // 4)), Acb(0.5, 0.6)]
            F_N_model = SRP.F_N_model(N₀, z)

            for inv_N in range(0, Arb(1 // N₀), 10)
                F_N_1 = F_N_model(inv_N)
                F_N_2 = SRP.hypgeom2f1(Acb(2inv_N), Acb(inv_N), Acb(1 + inv_N), z)

                @test Arblib.overlaps(F_N_1, F_N_2)
            end
        end
    end
end
