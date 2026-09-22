@testset "lemma 2.9" begin
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

    @testset "K" begin
        # Test that evaluation of the Taylor model agrees with direct
        # evaluation.
        N₀ = 64
        inv_N₀ = Arb(1 // N₀)
        z = SRP.exppii(Arb(1 // 4))
        t = z / 2

        res_1 = SRP.K_model(N₀, z)[2](t)(inv_N₀)

        res_2 = let
            ρ = SRP.c_N(inv_N₀)^2 * SRP.λ_app(inv_N₀)
            F_N_z = SRP.F_N_model(N₀, z)(inv_N₀)
            F_N_conj_z = SRP.F_N_model(N₀, conj(z))(inv_N₀)
            F_N_t = SRP.F_N_model(N₀, t)(inv_N₀)

            besselj0(
                sqrt(ρ) *
                abs(z)^inv_N₀ *
                sqrt(F_N_conj_z * (F_N_z - (t / z)^inv_N₀ * F_N_t)),
            )
        end

        @test Arblib.overlaps(res_1, res_2)

        # Test that the terms in the Taylor model expansion agrees
        # with the d's
        res = SRP.K_model(N₀, z)[2](t)

        @test Arblib.overlaps(SRP.polynomial(res)[0], SRP.d(0, z, t))
        @test Arblib.overlaps(SRP.polynomial(res)[1], SRP.d(1, z, t))
        @test Arblib.overlaps(SRP.polynomial(res)[2], SRP.d(2, z, t))
        @test Arblib.overlaps(SRP.polynomial(res)[3], SRP.d(3, z, t))
    end

    @testset "integral_K_4_V_l" begin
        z = SRP.exppii(Arb(1 // 4))

        res_1 = SRP.integral_K_4_V_l(64, 1, z)
        # Slightly more precise
        res_2 = SRP.integral_K_4_V_l(128, 1, z, a = Arb(1e-6), atol = 1e-4)
        # Even more precise
        res_3 = SRP.integral_K_4_V_l(256, 1, z, a = Arb(1e-8), atol = 1e-6)

        # This will throw an error if they do not all intersect
        @test isfinite(Arblib.intersection(real(res_1), real(res_2), real(res_3)))
        @test isfinite(Arblib.intersection(imag(res_1), imag(res_2), imag(res_3)))

        for l = 2:4
            res_1 = SRP.integral_K_4_V_l(64, l, z)
            # Slightly more precise
            res_2 = SRP.integral_K_4_V_l(128, l, z, a = Arb(0.025), atol = 1.0)
            # Even more precise
            res_3 = SRP.integral_K_4_V_l(256, l, z, a = Arb(0.0025), atol = 0.1)

            # This will throw an error if they do not all intersect
            @test isfinite(Arblib.intersection(real(res_1), real(res_2), real(res_3)))
            @test isfinite(Arblib.intersection(imag(res_1), imag(res_2), imag(res_3)))
        end

        # Test with wide input. In this case we only test for l = 1
        # since testing more cases is very time consuming.
        for t₀ in [1 / 3, 1 / 2, 0.99999, 1]
            t = setball(Arb, t₀, 1e-3)
            z = SRP.exppii(t)
            zₗ = SRP.exppii(Arblib.lbound(Arb, t))
            zᵤ = SRP.exppii(Arblib.ubound(Arb, t))

            # Use a higher N₀ so that we can get better enclosures
            # that increase the chance that we catch any errors.
            res_1 = SRP.integral_K_4_V_l(1000, 1, z, a = Arb(1e-8), atol = 1e-6)
            res_2 = SRP.integral_K_4_V_l(1000, 1, zₗ, a = Arb(1e-8), atol = 1e-6)
            res_3 = SRP.integral_K_4_V_l(1000, 1, zᵤ, a = Arb(1e-8), atol = 1e-6)

            @test Arblib.overlaps(res_1, res_2)
            @test Arblib.overlaps(res_1, res_3)
        end
    end
end
