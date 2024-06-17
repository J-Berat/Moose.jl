# Unit Tests
@testset "HighPass tests" begin
    # Test case 1: Simple known image and kernel
    image = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    kernel = ones(3, 3) / 9
    expected_output = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]  # Known result for this test case
    @test HighPass(image, kernel) ≈ expected_output

    # Test case 2: Random image and kernel
    image = rand(256, 256)
    kernel = ones(3, 3) / 9
    high_passed_image = HighPass(image, kernel)
    @test size(high_passed_image) == size(image)
    @test typeof(high_passed_image) == typeof(image)

    # Test case 3: Check if high-frequency components are enhanced
    image = rand(256, 256)
    kernel = ones(5, 5) / 25
    high_passed_image = HighPass(image, kernel)
    @test maximum(high_passed_image) > maximum(imfilter(image, kernel))

    # Test case 4: Edge case with zero image
    image = zeros(256, 256)
    kernel = ones(3, 3) / 9
    high_passed_image = HighPass(image, kernel)
    @test high_passed_image == zeros(256, 256)

    # Test case 5: Edge case with a kernel that is an impulse filter
    image = rand(256, 256)
    kernel = [0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]
    high_passed_image = HighPass(image, kernel)
    @test high_passed_image ≈ zeros(256, 256)
end
