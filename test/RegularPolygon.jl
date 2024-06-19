@testset "RegularPolygon" begin
    for N = 3:12
        domain = SRP.RegularPolygon(N)

        for i = 1:N
            @test SRP.polar_vertex(domain, SRP.vertex(domain, i + 1), i).φ ≈ 0 atol = 1e-15

            @test SRP.polar_vertex(domain, SRP.vertex(domain, i - 1), i).φ ≈
                  SRP.angle(domain, i)
        end
    end
end
