using Test
using Logging
using MOOSE

@testset "RMS" begin
    x = [1.0, 2.0, 3.0]
    @test isapprox(MOOSE.RMS(x), sqrt(2 / 3))

    x2 = [1.0, 2.0]
    y2 = [3.0, 4.0]
    @test isapprox(MOOSE.RMS(x2, y2), sqrt(0.5))

    z2 = [5.0, 6.0]
    @test isapprox(MOOSE.RMS(x2, y2, z2), sqrt(0.75))
end

@testset "Pnu" begin
    q = [1.0, 2.0]
    u = [3.0, 4.0]
    expected = sqrt.(q .^ 2 .+ u .^ 2)
    @test MOOSE.Pnu(q, u) == expected
end

@testset "Logging" begin
    Q = [0.1, 0.2]
    U = [0.3, 0.4]
    nu = [1.0, 1.1]
    phi = [-1.0, 0.0]

    buf = IOBuffer()
    logger = SimpleLogger(buf, Logging.Info)

    Logging.with_logger(logger) do
        MOOSE.RMSynthesis(Q, U, nu, phi; log_progress = true)
    end

    logs = String(take!(buf))
    @test occursin("RM synthesis complete", logs)
end
