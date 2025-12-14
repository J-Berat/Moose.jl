using Test
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

@testset "Input validation" begin
    tmp = mktempdir()
    missing_dir = joinpath(tmp, "missing")
    @test occursin("does not exist", MOOSE.ensure_directory_access(missing_dir))

    file_path = joinpath(tmp, "example.txt")
    open(file_path, "w") do io
        write(io, "sample")
    end

    msg = MOOSE.ensure_readable_file(file_path; expected_exts=[".fits"])
    @test occursin("expected extension", msg)

    nested_dir = joinpath(tmp, "folder")
    mkpath(nested_dir)
    msg_dir = MOOSE.ensure_readable_file(nested_dir; expected_exts=[".fits"])
    @test occursin("not a regular file", msg_dir)
end

@testset "ReadSimulation validation" begin
    err = nothing
    try
        MOOSE.ReadSimulation("/path/that/does/not/exist", "z", 1.0, 1.0, 1.0, 1.0)
    catch e
        err = e
    end
    @test err !== nothing
    @test occursin("does not exist", sprint(showerror, err))

    mktempdir() do dir
        for name in ["Bx", "By", "Bz", "density", "temperature", "Vx", "Vy", "Vz"]
            touch(joinpath(dir, "$(name).fits"))
        end

        corrupted_err = nothing
        try
            MOOSE.ReadSimulation(dir, "z", 1.0, 1.0, 1.0, 1.0)
        catch e
            corrupted_err = e
        end

        @test corrupted_err !== nothing
        @test occursin("could not be read", sprint(showerror, corrupted_err))
    end
end
