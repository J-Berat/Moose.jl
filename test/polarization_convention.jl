@testset "Ordered sky components and electric-vector angle" begin
    @test Moose.los_basis(:Ax, :Ay, :Az, "z") == (:Ax, :Ay, :Az)
    @test Moose.los_basis(:Ax, :Ay, :Az, "x") == (:Ay, :Az, :Ax)
    @test Moose.los_basis(:Ax, :Ay, :Az, "y") == (:Ax, :Az, :Ay)
    @test_throws ErrorException Moose.los_basis(1, 2, 3, "w")
    for T in (Float32, Float64)
        # Reference-axis B, second-axis B, and ±45-degree B orientations.
        b1 = T[1, 0, 1, 1, -1]
        b2 = T[0, 1, 1, -1, -1]
        psi = Moose.IntrinsicAngle(b1, b2)
        @test eltype(psi) == T
        @test cos.(2 .* psi) ≈ T[-1, 1, 0, 0, 0] atol=20eps(T)
        @test sin.(2 .* psi) ≈ T[0, 0, -1, 1, -1] atol=20eps(T)
        # The electric vector is perpendicular to B, including oblique fields.
        @test all(abs.(b1 .* cos.(psi) .+ b2 .* sin.(psi)) .< 20eps(T))
    end
    # Physical Cartesian components, passed through the same basis as both
    # full-cube and tiled pipelines. Expected angles are expressed directly
    # in physical components, independently of the helper's argument names.
    bx, by, bz = [2.0], [3.0], [5.0]
    for (los, expected) in (("x", atan(5.0, 3.0) + pi/2),
                            ("y", atan(5.0, 2.0) + pi/2),
                            ("z", atan(3.0, 2.0) + pi/2))
        b1, b2, _ = Moose.los_basis(bx, by, bz, los)
        @test only(Moose.IntrinsicAngle(b1, b2)) ≈ expected
    end
end
