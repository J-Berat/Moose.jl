@testset "max tests" begin
    cube = rand(100, 100, 50)
    max_values = max(cube)
    @test size(max_values) == (100, 100)
    @test maximum(max_values) == maximum(cube)
end

@testset "intLOS tests" begin
    cube = rand(100, 100, 50)
    PixelLength_cm = 1.0
    integrated_los = intLOS(cube, PixelLength_cm)
    @test size(integrated_los) == (100, 100)
    @test integrated_los[1, 1] == sum(cube[1, 1, :]) * PixelLength_cm
end

@testset "MeanSpectrum tests" begin
    cube = rand(100, 100, 50)
    mean_spectrum = MeanSpectrum(cube)
    @test size(mean_spectrum) == (50,)
    @test mean_spectrum[1] == mean(cube[:, :, 1])
end

@testset "MaxIndicesMap tests" begin
    cube = rand(100, 100, 50)
    ValueArray = range(-500, stop=500, length=50)
    MapMaxIndices = MaxIndicesMap(cube, ValueArray)
    @test size(MapMaxIndices) == (100, 100)
    max_indices = dropdims(argmax(cube, dims=3), dims=3)
    @test all(MapMaxIndices .== ValueArray[max_indices])
end

# The ask_user function cannot be tested via standard unit tests since it requires user interaction.

@testset "contains_fits_files tests" begin
    # Setup a temporary directory for testing
    temp_dir = mktempdir()
    try
        open(joinpath(temp_dir, "file.fits"), "w") do f end
        @test contains_fits_files(temp_dir) == true
        rm(joinpath(temp_dir, "file.fits"))
        @test contains_fits_files(temp_dir) == false
    finally
        rm(temp_dir, force=true)
    end
end

@testset "get_simulation_list tests" begin
    # Setup temporary directories for testing
    temp_base_dir = mktempdir()
    try
        dir1 = joinpath(temp_base_dir, "simu1")
        dir2 = joinpath(temp_base_dir, "simu2")
        mkpath(dir1)
        mkpath(dir2)
        open(joinpath(dir1, "file.fits"), "w") do f end
        simulation_list = get_simulation_list(temp_base_dir)
        @test length(simulation_list) == 1
        @test simulation_list[1] == dir1
    finally
        rm(temp_base_dir, force=true)
    end
end

@testset "display_simulations tests" begin
    simu_list = ["simu1", "simu2", "simu3"]
    # Capture the output of display_simulations
    io = IOBuffer()
    redirect_stdout(io) do
        display_simulations(simu_list)
    end
    output = String(take!(io))
    expected_output = "Available simulations:\n[1] simu1\n[2] simu2\n[3] simu3\n"
    @test output == expected_output
end


@testset "print_progress tests" begin
    io = IOBuffer()
    redirect_stdout(io) do
        print_progress(25, 100)
    end
    output = String(take!(io))
    expected_output = "\rProgress: |██████████                                        | 25/100"
    @test output == expected_output
end