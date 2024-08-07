@testset "RegularPolygon" begin
    @testset "polar_vertex" begin
        for N = 3:12
            domain = SRP.RegularPolygon{Arb}(N)

            for i = 1:N
                # Compute coordinates for all other vertices
                ps = map(i .+ (1:domain.N-1)) do j
                    SRP.polar_vertex(domain, SRP.vertex(domain, j), i)
                end

                # Angles should start at zero, end at the angle of the
                # vertex and be increasing
                @test Arblib.contains_zero(ps[1].φ)
                @test Arblib.overlaps(ps[end].φ, SRP.angle(domain, i))
                @test all(j -> ps[j].φ < ps[j+1].φ, 1:lastindex(ps)-1)

                # Radius should be increasing for first half and mirrored
                for j = 1:lastindex(ps)÷2-1
                    @test ps[j].r < ps[j+1].r
                end
                for j = 1:lastindex(ps)÷2
                    @test Arblib.overlaps(ps[j].r, ps[end-(j-1)].r)
                end

                # Compute random interior points
                ps = map(SRP.interior_points_random(domain, 20)) do xy
                    SRP.polar_vertex(domain, xy, i)
                end

                # Angles should be positive and less than the angle of
                # the vertex. The radius should be positive and less
                # than twice the radius of the vertex.
                for p in ps
                    @test 0 < p.φ < SRP.angle(domain, i)
                    @test 0 < p.r < 2SRP.Polar(SRP.vertex(domain, i)).r
                end
            end
        end
    end
end
