@testset "basic_functions" begin
    @testset "exppii" begin
        @test Arblib.overlaps(SRP.exppii(Arb(0.5)), SRP.exppii(Acb(0.5)))
    end

    @testset "abspow" begin
        for x0 in range(-2, 2, 21)
            for rad in logrange(1e-10, 4, 10)
                x = Arblib.setball(Arb, x0, rad)
                x_l, x_u = getinterval(Arb, x)
                for y in range(Arb(-2), 2, 21)
                    @test Arblib.overlaps(SRP.abspow(x, y), abs(x_l)^y)
                    @test Arblib.overlaps(SRP.abspow(x, y), abs(x_u)^y)
                end
            end
        end
    end

    @testset "logabspow" begin
        for x0 in range(-2, 2, 21)
            for rad in logrange(1e-10, 4, 10)
                x = Arblib.setball(Arb, x0, rad)
                x_l, x_u = getinterval(Arb, x)
                for y in range(Arb(-2), 2, 21)
                    for m = -1:4
                        @test Arblib.overlaps(
                            SRP.logabspow(x, m, y),
                            log(abs(x_l))^m * abs(x_l)^y,
                        )
                        @test Arblib.overlaps(
                            SRP.logabspow(x, m, y),
                            log(abs(x_u))^m * abs(x_u)^y,
                        )
                    end
                end
            end
        end
    end

    @testset "logpow" begin
        for re_0 in range(-2, 2, 21)
            for im_0 in range(-2, 2, 21)
                for rad in logrange(1e-10, 4, 10)
                    z = Acb(Arblib.setball(Arb, re_0, rad), Arblib.setball(Arb, im_0, rad))
                    re_l, re_u = getinterval(Arb, real(z))
                    im_l, im_u = getinterval(Arb, imag(z))
                    z_ll = Acb(re_l, im_l)
                    z_lu = Acb(re_l, im_u)
                    z_ul = Acb(re_u, im_l)
                    z_uu = Acb(re_u, im_u)
                    for y in range(Arb(-2), 2, 21)
                        for m = 0:4
                            @test Arblib.overlaps(SRP.logpow(z, m, y), log(z_ll)^m * z_ll^y)
                            @test Arblib.overlaps(SRP.logpow(z, m, y), log(z_lu)^m * z_lu^y)
                            @test Arblib.overlaps(SRP.logpow(z, m, y), log(z_ul)^m * z_ul^y)
                            @test Arblib.overlaps(SRP.logpow(z, m, y), log(z_uu)^m * z_uu^y)
                        end
                    end
                end
            end
        end
    end
end
