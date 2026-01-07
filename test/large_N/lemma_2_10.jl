@testset "lemma 2.10" begin
    @testset "integral_d_k_V_l" begin
        # Compute with different values for a and different tolerances and
        # verify that all results agree

        z = SRP.exppii(Arb(1 // 4))

        for (k, l) in [(1, 4), (2, 3), (3, 2)]
            res_1 = SRP.integral_d_k_V_l(k, l, z)
            # Slightly more precise
            res_2 = SRP.integral_d_k_V_l(k, l, z, a = Arb(1e-4), atol = 1e-5)
            # Even more precise
            res_3 = SRP.integral_d_k_V_l(k, l, z, a = Arb(1e-5), atol = 1e-6)

            # This will throw an error if they do not all intersect
            @test isfinite(Arblib.intersection(real(res_1), real(res_2), real(res_3)))
            @test isfinite(Arblib.intersection(imag(res_1), imag(res_2), imag(res_3)))
        end

        for (k, l) in [(2, 4), (3, 3), (3, 4)]
            res_1 = SRP.integral_d_k_V_l(k, l, z)
            # Slightly more precise
            res_2 = SRP.integral_d_k_V_l(k, l, z, a = Arb(0.05), atol = 1.0)
            # Even more precise
            res_3 = SRP.integral_d_k_V_l(k, l, z, a = Arb(0.005), atol = 0.1)

            # This will throw an error if they do not all intersect
            @test isfinite(Arblib.intersection(real(res_1), real(res_2), real(res_3)))
            @test isfinite(Arblib.intersection(imag(res_1), imag(res_2), imag(res_3)))
        end

        t = Arb((0, 1e-10))
        z = SRP.exppii(t)
        zᵤ = SRP.exppii(Arblib.ubound(Arb, t))

        for (k, l) in [(1, 4), (2, 3), (3, 2)]
            res_1 =
                SRP.integral_d_k_V_l(k, l, z, a = Arb(1e-5), b = Arb(0.9999), atol = 1e-6)
            res_2 = SRP.integral_d_k_V_l(k, l, zᵤ, a = Arb(1e-5), atol = 1e-6)

            @test Arblib.overlaps(res_1, res_2)
        end

        for (k, l) in [(2, 4), (3, 3), (3, 4)]
            res_1 =
                SRP.integral_d_k_V_l(k, l, z, a = Arb(0.005), b = Arb(0.999), atol = 0.1)
            res_2 = SRP.integral_d_k_V_l(k, l, zᵤ, a = Arb(0.005), atol = 0.1)

            @test Arblib.overlaps(res_1, res_2)
        end

        # Test with wide input. In this case we only test for (k, l) =
        # (2, 3) since testing more cases is very time consuming.
        for t₀ in [1 / 3, 1 / 2, 0.99999, 1]
            t = setball(Arb, t₀, 1e-4)
            z = SRP.exppii(t)
            zₗ = SRP.exppii(Arblib.lbound(Arb, t))
            zᵤ = SRP.exppii(Arblib.ubound(Arb, t))

            res_1 = SRP.integral_d_k_V_l(2, 3, z, a = Arb(1e-10), atol = 1e-7)
            res_2 = SRP.integral_d_k_V_l(2, 3, zₗ, a = Arb(1e-10), atol = 1e-8)
            res_3 = SRP.integral_d_k_V_l(2, 3, zᵤ, a = Arb(1e-10), atol = 1e-8)

            @test Arblib.overlaps(res_1, res_2)
            @test Arblib.overlaps(res_1, res_3)
        end
    end
end
