using Test
using MOOSE

@testset "CreateFreqFile" begin
    mktempdir() do dir
        start_freq = 100.0
        end_freq = 200.0
        num_freq = 10

        MOOSE.CreateFreqFile(start_freq, end_freq, num_freq, dir)

        path = joinpath(dir, "FreqHz.txt")
        freqs = readlines(path)

        @test length(freqs) == num_freq

        parsed_freqs = parse.(Float64, freqs)
        @test parsed_freqs[end] ≈ end_freq
    end
end
