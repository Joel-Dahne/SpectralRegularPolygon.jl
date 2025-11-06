@testset "basic_integrals" begin
    @testset "integral_log and integral_log_z" begin
        # Check that the difference when integrating to a and when
        # integrating to b is the same as the integral from a to b.

        a = Arb(0.25)
        b = Arb(0.5)
        z = Acb(0.4, 0.6)
        for m = 0:5
            @test Arblib.overlaps(
                SRP.integral_log(m, b) - SRP.integral_log(m, a),
                Arb(Arblib.integrate(x -> log(x)^m, a, b)),
            )

            @test Arblib.overlaps(
                SRP.integral_log_z(m, z, b) - SRP.integral_log_z(m, z, a),
                Arblib.integrate(t -> log(t / z)^m, a * z, b * z),
            )
        end
    end

    @testset "integral_log_1mtz and integral_logpow_1mtz" begin
        # Check that the difference when integrating from a and when
        # integrating from b is the same as the integral from a to b.

        a = Arb(0.25)
        b = Arb(0.5)
        z = Acb(0.4, 0.6)
        y = Arb(1.5)
        for m = 0:5
            @test Arblib.overlaps(
                SRP.integral_log_1mtz(z, m, a) - SRP.integral_log_1mtz(z, m, b),
                Arblib.integrate(t -> log(1 - t * z)^m, a, b),
            )

            @test Arblib.overlaps(
                SRP.integral_logpow_1mtz(z, m, y, a) - SRP.integral_logpow_1mtz(z, m, y, b),
                Arblib.integrate(
                    (t; analytic) ->
                        log(1 - t * z)^m *
                        Arblib.pow_analytic!(zero(t), 1 - t * z, Acb(y), analytic),
                    a,
                    b,
                    check_analytic = true,
                ),
            )
        end
    end
end
