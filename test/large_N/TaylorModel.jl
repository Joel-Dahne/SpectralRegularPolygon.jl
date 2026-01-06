@testset "$TaylorModel" for TaylorModel in [SRP.ArbTaylorModel, SRP.AcbTaylorModel]
    begin
        TaylorModel = SRP.ArbTaylorModel

        fs = [sin, cos, exp, atan]

        overlaps(x::Arb, y::Arb) = Arblib.overlaps(x, y)
        overlaps(x::Acb, y::Acb) = Arblib.overlaps(x, y)
        overlaps(x, y) =
            Arblib.overlaps(real(x), real(y)) && Arblib.overlaps(imag(x), imag(y))

        # Test if M looks like it is a Taylor model of f. We check so that
        # M overlaps with f at the midpoint, the endpoints and two points
        # close to the midpoint.
        test_f = (M, f) -> begin
            overlaps(M.p[0], f(M.x0)) || return false
            if !Arblib.isexact(M.I)
                x_l, x_u = getinterval(Arb, M.I)
                overlaps(M(x_l), f(x_l)) || return false
                overlaps(M(x_u), f(x_u)) || return false
                # Points close to the midpoint
                x1 = (4x_l + 6x_u) / 10
                x2 = (6x_l + 4x_u) / 10
                overlaps(M(x1), f(x1)) || return false
                overlaps(M(x2), f(x2)) || return false
            end
            return true
        end

        @testset "Construction" begin
            for f in fs
                for x0 in Arb[-2, 0, 1]
                    for r in Mag[0, 1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for n = 0:3
                            M1 = TaylorModel(f, I, x0, degree = n, enclosure_degree = -1)
                            M2 = TaylorModel(f, I, x0, degree = n, enclosure_degree = 0)
                            M3 = TaylorModel(f, I, x0, degree = n, enclosure_degree = 1)
                            @test Arblib.overlaps(M1, M2)
                            @test Arblib.overlaps(M1, M3)
                            @test Arblib.overlaps(M2, M3)
                            @test test_f(M1, f)
                            @test test_f(M2, f)
                            @test test_f(M3, f)
                        end
                    end
                end
            end

            p = TaylorModel == SRP.ArbTaylorModel ? ArbSeries((1, 2)) : AcbSeries((1, 2))

            M = zero(TaylorModel(p, Arb((-1, 1)), Arb(0.5)))
            @test iszero(M)
            @test !isone(M)
            @test isequal(M.I, Arb((-1, 1)))
            @test isequal(M.x0, 0.5)

            M = one(TaylorModel(p, Arb((-1, 1)), Arb(0.5)))
            @test isone(M)
            @test !iszero(M)
            @test isequal(M.I, Arb((-1, 1)))
            @test isequal(M.x0, 0.5)
        end

        @testset "checkcompatible" begin
            @test SRP.checkcompatible(
                Bool,
                TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
                TaylorModel(cos, Arb((-1, 1)), Arb(0), degree = 3),
            )

            @test !SRP.checkcompatible(
                Bool,
                TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
                TaylorModel(cos, Arb((-1, 2)), Arb(0), degree = 3),
            )

            @test !SRP.checkcompatible(
                Bool,
                TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
                TaylorModel(cos, Arb((-1, 1)), Arb(1), degree = 3),
            )

            @test !SRP.checkcompatible(
                Bool,
                TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
                TaylorModel(cos, Arb((-1, 1)), Arb(0), degree = 4),
            )

            @test_throws ErrorException SRP.checkcompatible(
                TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
                TaylorModel(cos, Arb((-1, 2)), Arb(0), degree = 3),
            )

            @test_throws ErrorException SRP.checkcompatible(
                TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
                TaylorModel(cos, Arb((-1, 1)), Arb(1), degree = 3),
            )

            @test_throws ErrorException SRP.checkcompatible(
                TaylorModel(sin, Arb((-1, 1)), Arb(0), degree = 3),
                TaylorModel(cos, Arb((-1, 1)), Arb(0), degree = 4),
            )
        end

        @testset "truncate" begin
            for f in fs
                for x0 in Arb[-2, 0, 1]
                    for r in Mag[0, 1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for n = 0:3
                            M = TaylorModel(f, I, x0, degree = n)
                            for m = 0:n
                                M_truncated = SRP.truncate(M, degree = m)
                                @test test_f(M_truncated, f)
                            end
                        end
                    end
                end
            end
        end

        @testset "compose" begin
            for f in fs
                for g in fs
                    for x0 in Arb[-2, 0, 1]
                        for r in Mag[0, 1e-10, 1e-5, 1e0]
                            I = add_error(x0, r)
                            for n = 0:3
                                Mf = TaylorModel(f, I, x0, degree = n)
                                Mgf = TaylorModel(g ∘ f, I, x0, degree = n)

                                @test Arblib.overlaps(Mgf, SRP.compose(g, Mf))
                                @test test_f(SRP.compose(g, Mf), g ∘ f)
                            end
                        end
                    end
                end
            end
        end

        @testset "Arithmetic" begin
            for f in fs
                for g in fs
                    for x0 in Arb[-2, 0, 1]
                        for r in Mag[0, 1e-10, 1e-5, 1e0]
                            I = add_error(x0, r)
                            for n = 0:3
                                Mf = TaylorModel(f, I, x0, degree = n)
                                Mg = TaylorModel(g, I, x0, degree = n)

                                Mf_plus_g = TaylorModel(x -> f(x) + g(x), I, x0, degree = n)
                                Mf_minus_g =
                                    TaylorModel(x -> f(x) - g(x), I, x0, degree = n)
                                Mf_mul_g = TaylorModel(x -> f(x) * g(x), I, x0, degree = n)
                                Mneg_f = TaylorModel(x -> -f(x), I, x0, degree = n)

                                @test Arblib.overlaps(Mf + Mg, Mf_plus_g)
                                @test Arblib.overlaps(Mf - Mg, Mf_minus_g)
                                @test Arblib.overlaps(Mf * Mg, Mf_mul_g)
                                @test Arblib.overlaps(-Mf, Mneg_f)

                                @test test_f(Mf + Mg, x -> f(x) + g(x))
                                @test test_f(Mf - Mg, x -> f(x) - g(x))
                                @test test_f(Mf * Mg, x -> f(x) * g(x))
                                @test test_f(-Mf, x -> -f(x))

                                if !Arblib.contains_zero(Mg.p[0])
                                    Mf_div_g =
                                        TaylorModel(x -> f(x) / g(x), I, x0, degree = n)

                                    @test Arblib.overlaps(Mf / Mg, Mf_div_g)
                                    @test test_f(Mf / Mg, x -> f(x) / g(x))
                                end
                            end
                        end
                    end
                end
            end
        end

        @testset "Scalar arithmetic" begin
            for f in fs
                for x0 in Arb[-2, 0, 1]
                    for r in Mag[0, 1e-10, 1e-5, 1e0]
                        # There is a bug in Arblib.jl which makes this case throw an error
                        TaylorModel == SRP.AcbTaylorModel && iszero(x0) && continue

                        I = add_error(x0, r)
                        for n = 0:3
                            Mf = TaylorModel(f, I, x0, degree = n)

                            c = Arb(2)
                            Mf_add_c = TaylorModel(x -> f(x) + c, I, x0, degree = n)
                            Mf_sub_c = TaylorModel(x -> f(x) - c, I, x0, degree = n)
                            Mf_mul_c = TaylorModel(x -> f(x) * c, I, x0, degree = n)
                            Mf_div_c = TaylorModel(x -> f(x) / c, I, x0, degree = n)
                            c_add_Mf = TaylorModel(x -> c + f(x), I, x0, degree = n)
                            c_sub_Mf = TaylorModel(x -> c - f(x), I, x0, degree = n)
                            c_mul_Mf = TaylorModel(x -> c * f(x), I, x0, degree = n)
                            c_div_Mf = TaylorModel(x -> c / f(x), I, x0, degree = n)

                            @test Arblib.overlaps(Mf + c, Mf_add_c)
                            @test Arblib.overlaps(Mf - c, Mf_sub_c)
                            @test Arblib.overlaps(Mf * c, Mf_mul_c)
                            @test Arblib.overlaps(Mf / c, Mf_div_c)
                            @test Arblib.overlaps(c + Mf, Mf_add_c)
                            @test Arblib.overlaps(c - Mf, c_sub_Mf)
                            @test Arblib.overlaps(c * Mf, Mf_mul_c)
                            @test Arblib.overlaps(c / Mf, c_div_Mf)
                        end
                    end
                end
            end
        end

        @testset "Shift" begin
            for f in fs
                for x0 in Arb[-2, 0, 1]
                    for r in Mag[0, 1e-10, 1e-5, 1e0]
                        I = add_error(x0, r)
                        for n = 0:3
                            Mf = TaylorModel(f, I, x0, degree = n)
                            Mf_mul_x =
                                TaylorModel(x -> (x - x0) * f(x), I, x0, degree = n + 1)
                            Mf_mul_x2 =
                                TaylorModel(x -> (x - x0)^2 * f(x), I, x0, degree = n + 2)

                            @test Arblib.overlaps(Mf, Mf_mul_x << 1)
                            @test Arblib.overlaps(Mf, Mf_mul_x2 << 2)
                            @test Arblib.overlaps(Mf_mul_x, Mf_mul_x2 << 1)
                            @test Arblib.overlaps(Mf >> 1, Mf_mul_x)
                            @test Arblib.overlaps(Mf >> 2, Mf_mul_x2)
                            @test Arblib.overlaps(Mf_mul_x >> 1, Mf_mul_x2)
                        end
                    end
                end
            end
        end

        if TaylorModel == SRP.AcbTaylorModel
            @testset "abs" begin
                for f in [identity, sin, exp, y -> -y]
                    for x0 in Arb[0.1, 1.6, 2.7]
                        for r in Mag[0, 1e-10, 1e-5]
                            I = add_error(x0, r)
                            for n = 0:3
                                My = SRP.AcbTaylorModel(f, I, x0, degree = n)
                                M = abs(My)
                                @test test_f(M, abs∘f)
                            end
                        end
                    end
                end
            end
        end
    end
end
