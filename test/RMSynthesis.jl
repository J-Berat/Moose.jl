# Unit Tests
@testset "RMSynthesis tests" begin
    Q = [0.1, 0.2, 0.3]
    U = [0.4, 0.5, 0.6]
    nuArray = [1e9, 1.1e9, 1.2e9]
    PhiArray = [-100, 0, 100]
    
    absF, realF, imagF = RMSynthesis(Q, U, nuArray, PhiArray)
    
    @test length(absF) == length(PhiArray)
    @test length(realF) == length(PhiArray)
    @test length(imagF) == length(PhiArray)
    
    @test all(absF .>= 0)
end

@testset "getRMSF tests" begin
    nuArray = [1e9, 1.1e9, 1.2e9]
    PhiArray = [-100, 0, 100]
    
    absRMSF, fwhmRMSF = getRMSF(nuArray, PhiArray)
    
    @test length(absRMSF) == length(PhiArray)
    @test fwhmRMSF > 0
    
    @test all(absRMSF .>= 0)
end