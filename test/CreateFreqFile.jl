@testset "CreateFreqFile tests" begin
    # Test case 1: Normal case
    CreateFreqFile(10.0, 100.0, 10, "test_dir")
    @test isfile("test_dir/FreqHz.txt")
    rm("test_dir/FreqHz.txt")

    # Test case 2: Negative frequency should throw an error
    @test_throws ArgumentError create_freq_file(-10.0, 100.0, 10, "test_dir")

    # Test case 3: Zero number of frequencies should throw an error
    @test_throws ArgumentError create_freq_file(10.0, 100.0, 0, "test_dir")

    # Test case 4: Empty repertory should create file in current directory
    create_freq_file(10.0, 100.0, 10)
    @test isfile("FreqHz.txt")
    rm("FreqHz.txt")
end