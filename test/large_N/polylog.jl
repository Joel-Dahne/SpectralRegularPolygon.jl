@testset "polylog" begin
    @testset "S" begin
        # z values we test for
        zs = [
            Acb(1 // 4),
            Acb(1 // 2, 1 // 3),
            SRP.exppii(Arb(1 // 4)),
            SRP.exppii(Arb(1 // 2)),
            SRP.exppii(Arb(3 // 4)),
        ]

        # For n = 2, S(2, z) is explicitly given by polylog(2, z) and
        # we test that they agree.
        for z in zs
            @test Arblib.overlaps(SRP.S(2, z), 2SRP.polylog(2, z))
        end

        # Verify that the compute values for S(n, z) agrees with some
        # values computed with Mathemtica. The Mathematica values have
        # been computed with
        # S3[z_] := 2^2*PolyLog[1, 2, z] - 2*PolyLog[2, 1, z]
        # S4[z_] := 2^3*PolyLog[1, 3, z] - 2^2*PolyLog[2, 2, z] + 2*PolyLog[3, 1, z]
        # S5[z_] := 2^4*PolyLog[1, 4, z] - 2^3*PolyLog[2, 3, z] + 2^2*PolyLog[3, 2, z] - 2*PolyLog[4, 1, z]
        # We compute 20 digits in Mathematica for all cases except
        # when z = im, when we compute 30 digits (since our enclosures
        # are more accurate for this case)
        S_values = [
            [
                Acb("-0.44181297095909025719"),
                Acb("-0.99607508949881116070", "-0.25856838466058147081"),
                Acb("-2.2228496984364014122", "-0.6753691477336351257"),
                Acb("-0.47347697528897258097810384350", "-2.39248796291835333473342769958"),
                Acb("1.5811903334693655311", "-1.8302124476431511538"),
            ],
            [
                Acb("0.48246890942873742499"),
                Acb("0.85266401101240961663", "0.61503955866660442686"),
                Acb("1.01103206541224357038", "0.67195702729321656807"),
                Acb("0.56187550449994927762384880517", "1.96771842432445228287908411985"),
                Acb("-1.4258490754270000146", "1.8769725893707754301"),
            ],
            [
                Acb("-0.48903357274388432581"),
                Acb("-1.00203338847885385577", "-0.63074492031658443153"),
                Acb("-1.2053325647507299970", "-1.4956868529609937333"),
                Acb("-0.22469407094175946035686851392", "-1.86516746284522705541617284324"),
                Acb("1.3650692803254877513", "-1.6783532510585516166"),
            ],
        ]

        for (n, ress) in zip(3:5, S_values)
            for (z, res) in zip(zs, ress)
                @test Arblib.contains(SRP.S(n, z), res)
            end
        end
    end

    @testset "Multiple polylogs" begin
        # We compare the values for the multiple polylog
        # implementations with values computed using a Mathematica
        # implementation. The implementation in Mathematica uses the
        # same expressions for the polylogarithms, but since it was
        # implemented separately this reduces the risk of accidental
        # typos.

        # z values we test for

        # TODO: Actually add some tests here!

        zs = []


        polylog_1_1_values = []
        polylog_1_2_values = []
        polylog_1_3_values = []
        polylog_2_1_values = []
        polylog_2_2_values = []
        polylog_3_1_values = []
        polylog_1_1_1_values = []
        polylog_1_1_2_values = []
        polylog_1_2_1_values = []
        polylog_2_1_1_values = []
        polylog_1_1_1_1_values = []

        for (z, res) in zip(zs, polylog_1_1_values)
            @test Arblib.contains(SRP.polylog_1_1(z), res)
        end
        for (z, res) in zip(zs, polylog_1_2_values)
            @test Arblib.contains(SRP.polylog_1_2(z), res)
        end
        for (z, res) in zip(zs, polylog_1_3_values)
            @test Arblib.contains(SRP.polylog_1_3(z), res)
        end
        for (z, res) in zip(zs, polylog_2_1_values)
            @test Arblib.contains(SRP.polylog_2_1(z), res)
        end
        for (z, res) in zip(zs, polylog_2_2_values)
            @test Arblib.contains(SRP.polylog_2_2(z), res)
        end
        for (z, res) in zip(zs, polylog_3_1_values)
            @test Arblib.contains(SRP.polylog_3_1(z), res)
        end
        for (z, res) in zip(zs, polylog_1_1_1_values)
            @test Arblib.contains(SRP.polylog_1_1_1(z), res)
        end
        for (z, res) in zip(zs, polylog_1_1_2_values)
            @test Arblib.contains(SRP.polylog_1_1_2(z), res)
        end
        for (z, res) in zip(zs, polylog_1_2_1_values)
            @test Arblib.contains(SRP.polylog_1_2_1(z), res)
        end
        for (z, res) in zip(zs, polylog_2_1_1_values)
            @test Arblib.contains(SRP.polylog_2_1_1(z), res)
        end
        for (z, res) in zip(zs, polylog_1_1_1_1_values)
            @test Arblib.contains(SRP.polylog_1_1_1_1(z), res)
        end
    end
end
