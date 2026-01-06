@testset "V" begin
    for re_0 in range(-0.6, 0.6, 11)
        for im_0 in range(-0.6, 0.6, 11)
            z = Acb(re_0, im_0)
            Ds = [SRP.V_div_z_bound(l, Arblib.abs_ubound(Arb, z)) for l = 1:4]
            for l = 1:4
                if abs(SRP.V(l, z) / z) > Ds[l]
                    @show re_0 im_0 l
                end

                @test !(abs(SRP.V(l, z) / z) > Ds[l])
            end
        end
    end
end
