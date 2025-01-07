"""
    ProcessSynchrotron(simu::String, LOS::String, FaradayRotation::String, responseSynchrotron::String, 
                       df::DataFrame, kernel_size_synchrotron::Float64, zeta::Float64, Geff::Float64, 
                       omegaPAH::Float64, XC::Float64, nuArray::AbstractArray, PhiArray::AbstractArray, 
                       PixelLength_pc::Float64, PixelLength_cm::Float64, BoxLength_pc::Float64, 
                       DistanceArray::AbstractArray)

Process synchrotron emission data for a given simulation and line of sight (LOS).

# Arguments
- `simu::String`: The base directory of the simulation.
- `LOS::String`: The line of sight direction ("x", "y", or "z").
- `FaradayRotation::String`: Flag to indicate whether to include Faraday rotation ("Y" or "N").
- `responseSynchrotron::String`: Flag to indicate whether to apply synchrotron response filtering ("Y" or "N").
- `df::DataFrame`: DataFrame containing additional data for synchrotron computation.
- `kernel_size_synchrotron::Float64`: Size of the Gaussian kernel for filtering.
- `zeta::Float64`: Parameter for electron density computation.
- `Geff::Float64`: Parameter for electron density computation.
- `omegaPAH::Float64`: Parameter for electron density computation.
- `XC::Float64`: Parameter for electron density computation.
- `nuArray::AbstractArray`: Array of frequency values.
- `PhiArray::AbstractArray`: Array of Faraday depth values.
- `PixelLength_pc::Float64`: Pixel length in parsecs.
- `PixelLength_cm::Float64`: Pixel length in centimeters.
- `BoxLength_pc::Float64`: Box length in parsecs.
- `DistanceArray::AbstractArray`: Array of distance values.

# Description
This function processes synchrotron emission data for a given simulation and line of sight (LOS). It reads simulation parameters, computes electron density, and calculates various quantities such as the integral of quantities, Faraday rotation, Stokes parameters Q and U, brightness temperature, synchrotron power, and performs RM synthesis if required. The results are saved in FITS files in the appropriate directories.

# Example
```julia
# Example usage
simu = "path/to/simulation"
LOS = "z"
FaradayRotation = "Y"
responseSynchrotron = "Y"
df = DataFrame()
kernel_size_synchrotron = 1.0
zeta = 1.0
Geff = 1.0
omegaPAH = 1.0
XC = 1.0
nuArray = [1e9, 2e9, 3e9]
PhiArray = [-100, 0, 100]
PixelLength_pc = 0.1
PixelLength_cm = 3.086e18
BoxLength_pc = 10.0
DistanceArray = [1.0, 2.0, 3.0]

ProcessSynchrotron(simu, LOS, FaradayRotation, responseSynchrotron, df, kernel_size_synchrotron, zeta, Geff, omegaPAH, XC, nuArray, PhiArray, PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray)
"""
function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String, 
                       df::DataFrame, kernel_size_synchrotron, zeta::Float64, Geff::Float64, 
                       omegaPAH::Float64, XC::Float64, nuArray::AbstractArray, PhiArray, 
                       PixelLength_pc::Float64, PixelLength_cm::Float64, BoxLength_pc, 
                       DistanceArray::AbstractArray, conversionn, conversionT, conversionB)
    
    #default_path = joinpath(simu, LOS, "Synchrotron", "zeta-16")
    #resultspath = ask_user("Where do you want to save your files?", default_path)
    #mkpath(resultspath) 
    println("Processing Synchrotron data for LOS: $LOS")
    resultspath = joinpath(simu, LOS, "Synchrotron")
    mkpath(resultspath)  
    # Read simulation parameters
    SimuParameters = ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    B1, B2, BLOS = SimuParameters[1], SimuParameters[2], SimuParameters[3]
    T, n = SimuParameters[4], SimuParameters[5]
    SimuParameters = nothing
    
    Bperpcube = Bperp(B1, B2)
    psi_src = IntrinsicAngle(B1, B2)

    # Compute electron density
    println("Computing electron density")
    ne = Wolfire_ne(zeta, Geff, omegaPAH, XC, T, n)
    WriteData3D(resultspath, ne, "ne", DistanceArray)

    # Compute integral of quantities
    println("Computing integrated quantities")
    Btotal = Btot(B1,B2,BLOS)
    intBtotal = intLOS(Btotal, PixelLength_cm)
    sigmaBtotal = sigmaLOS(Btotal)
    Btotal = nothing
    Ne = intLOS(ne, PixelLength_cm)
    sigmane = sigmaLOS(ne)
    sigmaT = sigmaLOS(T)
    intBLOS = intLOS(BLOS, PixelLength_cm)
    sigmaBLOS = sigmaLOS(BLOS)
    intBperp = intLOS(Bperpcube, PixelLength_cm)
    B1 = nothing
    B2 = nothing
    
    WriteData2D(resultspath, intBtotal, "intBtotal")
    WriteData2D(resultspath, sigmaBtotal, "sigmaBtotal")
    WriteData2D(resultspath, Ne, "intne")
    WriteData2D(resultspath, sigmane, "sigmane")
    WriteData2D(resultspath, sigmaT, "sigmaT")
    WriteData2D(resultspath, intBLOS, "intBLOS")
    WriteData2D(resultspath, sigmaBLOS, "sigmaBLOS")
    WriteData2D(resultspath, intBperp, "intBperp")

    # Faraday rotation
    if uppercase(FaradayRotation) == "Y"
        resultspath = joinpath(resultspath, "WithFaraday")
        mkpath(resultspath)
        println("Computing RM")
        dRM = deltaRM(BLOS, ne, PixelLength_pc)
        RMcube = RM(dRM)
        RMmap = RMcube[:, :, end]
        BLOS = nothing
        WriteData2D(resultspath, RMmap, "RMmap")
    else
        resultspath = joinpath(resultspath, "noFaraday")
        mkpath(resultspath)
        println("No Faraday rotation included.")
    end

    # Compute Q and U
    if uppercase(FaradayRotation) == "Y"
        println("Computing Qnu and Unu with Faraday rotation")
        Qnu, Unu = QUnu3D(Bperpcube, psi_src, RMcube, nuArray, df, PixelLength_cm)
    else
        println("Computing Qnu and Unu without Faraday rotation")
        Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, PixelLength_cm)
    end
    T_nu = Tnu3D(Bperpcube, nuArray, df, PixelLength_cm)
    Bperpcube = nothing

    if uppercase(responseSynchrotron) == "Y"
        kernel = Kernel.gaussian(kernel_size_synchrotron)
        println("Applying filtering")
        for i in 1:length(nuArray)
            Qnu[:, :, i] = HighPass(Qnu[:, :, i], kernel)
            Unu[:, :, i] = HighPass(Unu[:, :, i], kernel)
            T_nu[:, :, i] = HighPass(T_nu[:, :, i], kernel)
        end
        resultspath = joinpath(resultspath, "filtered")
        mkpath(resultspath)
    else
        println("No filtering performed.")
    end

    WriteData3D(resultspath, Qnu, "Qnu", nuArray)
    WriteData3D(resultspath, Unu, "Unu", nuArray)
    WriteData3D(resultspath, T_nu, "T_nu", nuArray)
    mv(joinpath(resultspath, "T_nu.fits"), joinpath(resultspath, "Tnu.fits"),force=true)  # Rename file

    println("Computing Pnu and Pnumax")
    Pnucube = Pnu(Qnu, Unu)
    Pnumax = max(Pnucube)
    WriteData3D(resultspath, Pnucube, "Pnu", nuArray)
    WriteData2D(resultspath, Pnumax, "Pnumax")

    if uppercase(FaradayRotation) == "Y"
        println("Performing RMsynthesis")
        FDF, realFDF, imagFDF = RMSynthesis(Qnu, Unu, nuArray * 1e6, PhiArray)
        Pmax = max(FDF)
        WriteData3D(resultspath, FDF, "FDF", PhiArray)
        WriteData3D(resultspath, realFDF, "realFDF", PhiArray)
        WriteData3D(resultspath, imagFDF, "imagFDF", PhiArray)
        WriteData2D(resultspath, Pmax, "Pmax")
    else
        println("No Faraday tomography performed.")
    end

    #if uppercase(MomentsMaps) == "Y"
    #    println("Computing Faraday moments maps")
    #    M0, M1, M2 = moments_map(FDF, PhiArray)
    #    WriteData2D(resultspath, M0, "M0_Faraday")
    #    WriteData2D(resultspath, M1, "M1_Faraday")
    #    WriteData2D(resultspath, M2, "M2_Faraday")
    #else
    #    println("No Moment computed.")
    #end
end

function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String, 
                       df::DataFrame, kernel_size_synchrotron, nuArray::AbstractArray, PhiArray, 
                       PixelLength_pc::Float64, PixelLength_cm::Float64, BoxLength_pc, 
                       DistanceArray::AbstractArray, conversionn, conversionT, conversionB)
    
    println("Processing Synchrotron data for LOS: $LOS")
    resultspath = joinpath(simu, LOS, "Synchrotron")
    mkpath(resultspath)  # Create directory if it doesn't exist
    
    # Read simulation parameters
    SimuParameters = ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    B1, B2, BLOS = SimuParameters[1], SimuParameters[2], SimuParameters[3]
    T, n = SimuParameters[4], SimuParameters[5]
    nH2, nHp = SimuParameters[6], SimuParameters[7]
    SimuParameters = nothing
    
    Bperpcube = Bperp(B1, B2)
    psi_src = IntrinsicAngle(B1, B2)

    # Compute integral of quantities
    println("Computing integrated quantities")
    Btotal = Btot(B1,B2,BLOS)
    intBtotal = intLOS(Btotal, 1.0)
    sigmaBtotal = sigmaLOS(Btotal)
    Btotal = nothing
    Ne = intLOS(nHp, PixelLength_cm)
    sigmane = sigmaLOS(nHp)
    sigmaT = sigmaLOS(T)
    intBLOS = intLOS(BLOS, PixelLength_cm)
    sigmaBLOS = sigmaLOS(BLOS)
    intBperp = intLOS(Bperpcube, PixelLength_cm)
    B1 = nothing
    B2 = nothing
    
    WriteData2D(resultspath, intBtotal, "intBtotal")
    WriteData2D(resultspath, sigmaBtotal, "sigmaBtotal")
    WriteData2D(resultspath, Ne, "intne")
    WriteData2D(resultspath, sigmane, "sigmane")
    WriteData2D(resultspath, sigmaT, "sigmaT")
    WriteData2D(resultspath, intBLOS, "intBLOS")
    WriteData2D(resultspath, sigmaBLOS, "sigmaBLOS")
    WriteData2D(resultspath, intBperp, "intBperp")

    # Faraday rotation
    if uppercase(FaradayRotation) == "Y"
        resultspath = joinpath(resultspath, "WithFaraday")
        mkpath(resultspath)
        println("Computing RM")
        dRM = deltaRM(BLOS, nHp, PixelLength_pc)
        RMcube = RM(dRM)
        RMmap = RMcube[:, :, end]
        BLOS = nothing
        WriteData2D(resultspath, RMmap, "RMmap")
    else
        resultspath = joinpath(resultspath, "noFaraday")
        mkpath(resultspath)
        println("No Faraday rotation included.")
    end

    # Compute Q and U
    if uppercase(FaradayRotation) == "Y"
        println("Computing Qnu and Unu with Faraday rotation")
        Qnu, Unu = QUnu3D(Bperpcube, psi_src, RMcube, nuArray, df, PixelLength_cm)
    else
        println("Computing Qnu and Unu without Faraday rotation")
        Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, PixelLength_cm)
    end
    T_nu = Tnu3D(Bperpcube, nuArray, df, PixelLength_cm)
    Bperpcube = nothing

    if uppercase(responseSynchrotron) == "Y"
        kernel = Kernel.gaussian(kernel_size_synchrotron)
        println("Applying filtering")
        for i in 1:length(nuArray)
            Qnu[:, :, i] = HighPass(Qnu[:, :, i], kernel)
            Unu[:, :, i] = HighPass(Unu[:, :, i], kernel)
            T_nu[:, :, i] = HighPass(T_nu[:, :, i], kernel)
        end
        resultspath = joinpath(resultspath, "filtered")
        mkpath(resultspath)
    else
        println("No filtering performed.")
    end

    WriteData3D(resultspath, Qnu, "Qnu", nuArray)
    WriteData3D(resultspath, Unu, "Unu", nuArray)
    WriteData3D(resultspath, T_nu, "T_nu", nuArray)
    mv(joinpath(resultspath, "T_nu.fits"), joinpath(resultspath, "Tnu.fits"),force=true)  # Rename file

    println("Computing Pnu and Pnumax")
    Pnucube = Pnu(Qnu, Unu)
    Pnumax = max(Pnucube)
    WriteData3D(resultspath, Pnucube, "Pnu", nuArray)
    WriteData2D(resultspath, Pnumax, "Pnumax")

    if uppercase(FaradayRotation) == "Y"
        println("Performing RMsynthesis")
        FDF, realFDF, imagFDF = RMSynthesis(Qnu, Unu, nuArray * 1e6, PhiArray)
        Pmax = max(FDF)
        WriteData3D(resultspath, FDF, "FDF", PhiArray)
        WriteData3D(resultspath, realFDF, "realFDF", PhiArray)
        WriteData3D(resultspath, imagFDF, "imagFDF", PhiArray)
        WriteData2D(resultspath, Pmax, "Pmax")
    else
        println("No Faraday tomography performed.")
    end
end