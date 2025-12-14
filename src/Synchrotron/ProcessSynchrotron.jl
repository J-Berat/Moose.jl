using DataFrames

function compute_los_spacing(box_length_pc, depth)
    depth > 0 || error("Cannot compute line-of-sight spacing for an empty cube.")
    los_pixel_length_pc = Float64(box_length_pc) / depth
    los_pixel_length_cm = los_pixel_length_pc * PARSEC_TO_CM
    los_distance_array = range(start=0.0, step=los_pixel_length_pc, length=depth)
    return los_pixel_length_pc, los_pixel_length_cm, los_distance_array
end
function ProcessSynchrotron(simu::AbstractString, LOS, FaradayRotation::AbstractString, responseSynchrotron::AbstractString,
                       df::DataFrame, add_noise, Noise_nu, kernel_size_synchrotron, zeta::Float64, Geff::Float64,
                       omegaPAH::Float64, XC::Float64, nuArray::AbstractArray, PhiArray,
                       PixelLength_pc, PixelLength_cm, BoxLength_pc,
                       DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false)
    #default_path = joinpath(simu, LOS, "Synchrotron")
    #resultspath = ask_user("Where do you want to save your files?", default_path)
    #mkpath(resultspath) 

    #println("Processing Synchrotron data for LOS: $LOS")
    resultspath = joinpath(simu, LOS, "Synchrotron")
    mkpath(resultspath)  
    
    # Read simulation parameters
    SimuParameters = ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    B1, B2, BLOS = SimuParameters[1], SimuParameters[2], SimuParameters[3]
    T, n = SimuParameters[4], SimuParameters[5]
    SimuParameters = nothing

    los_pixel_length_pc, los_pixel_length_cm, los_distance_array = compute_los_spacing(BoxLength_pc, size(B1, 3))
    PixelLength_pc, PixelLength_cm, DistanceArray =
        los_pixel_length_pc, los_pixel_length_cm, los_distance_array

    Bperpcube = Bperp(B1, B2)
    psi_src = IntrinsicAngle(B1, B2)

    faraday_enabled = uppercase(FaradayRotation) == "Y"
    filtering_enabled = uppercase(responseSynchrotron) == "Y"
    noise_enabled = uppercase(add_noise) == "Y"

    # Compute electron density
    println("-------------------------------------------")
    println(Crayon(foreground=:blue, bold=true)("Computing electron density"))
    ne = Wolfire_ne(zeta, Geff, omegaPAH, XC, T, n)
    WriteData3D(resultspath, ne, "ne", los_distance_array)

    # Compute integral of quantities
    println("-------------------------------------------")
    println(Crayon(foreground=:blue, bold=true)("Computing integrated quantities"))
    Btotal = Btot(B1,B2,BLOS)
    intBtotal = intLOS(Btotal, los_pixel_length_cm)
    sigmaBtotal = sigmaLOS(Btotal)
    Btotal = nothing
    Ne = intLOS(ne, los_pixel_length_cm)
    sigmane = sigmaLOS(ne)
    sigmaT = sigmaLOS(T)
    intBLOS = intLOS(BLOS, los_pixel_length_cm)
    sigmaBLOS = sigmaLOS(BLOS)
    intBperp = intLOS(Bperpcube, los_pixel_length_cm)
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
    if faraday_enabled
        resultspath = joinpath(resultspath, "WithFaraday")
        mkpath(resultspath)
        println("-------------------------------------------")
        println(Crayon(foreground=:blue, bold=true)("Computing RM"))
        dRM = deltaRM(BLOS, ne, los_pixel_length_pc)
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
    if faraday_enabled
        println("-------------------------------------------")
        println(Crayon(foreground=:blue, bold=true)("Computing Qnu and Unu with Faraday rotation"))
        Qnu, Unu = QUnu3D(Bperpcube, psi_src, RMcube, nuArray, df, los_pixel_length_cm)
    else
        println("-------------------------------------------")
        println(Crayon(foreground=:blue, bold=true)("Computing Qnu and Unu without Faraday rotation"))
        Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, los_pixel_length_cm)
    end
    T_nu = Tnu3D(Bperpcube, nuArray, df, los_pixel_length_cm)
    Bperpcube = nothing
    
    if filtering_enabled
        kernel = Kernel.gaussian(kernel_size_synchrotron)
        println("-------------------------------------------")
        println(Crayon(foreground=:red, bold=true)("Applying filtering"))
        @views for i in axes(Qnu, 3)
            Qnu[:, :, i] = HighPass(Qnu[:, :, i], kernel)
            Unu[:, :, i] = HighPass(Unu[:, :, i], kernel)
            T_nu[:, :, i] = HighPass(T_nu[:, :, i], kernel)
        end
        resultspath = joinpath(resultspath, "filtered")
        mkpath(resultspath)
    else
        println("No filtering performed.")
    end
    # Adding noise
    if noise_enabled
        println("-------------------------------------------")
        println(Crayon(foreground=:red, bold=true)("Adding gaussian noise to Q and U with sigma = $Noise_nu"))
        noiseQ = similar(Qnu[:, :, 1])
        noiseU = similar(Unu[:, :, 1])
        distQ = Normal(0, Noise_nu)
        distU = Normal(0, Noise_nu)
        @views for i in axes(Qnu, 3)
            rand!(distQ, noiseQ)
            rand!(distU, noiseU)
            Qnu[:, :, i] .+= noiseQ
            Unu[:, :, i] .+= noiseU
        end
    end

    WriteData3D(resultspath, Qnu, "Qnu", nuArray)
    WriteData3D(resultspath, Unu, "Unu", nuArray)
    WriteData3D(resultspath, T_nu, "T_nu", nuArray)
    mv(joinpath(resultspath, "T_nu.fits"), joinpath(resultspath, "Tnu.fits"),force=true)  # Rename file
    
    println("-------------------------------------------")
    println(Crayon(foreground=:blue, bold=true)("Computing Pnu and Pnumax"))
    Pnucube = Pnu(Qnu, Unu)
    Pnumax = maxCube(Pnucube)
    WriteData3D(resultspath, Pnucube, "Pnu", nuArray)
    WriteData2D(resultspath, Pnumax, "Pnumax")

    if faraday_enabled
        println("-------------------------------------------")
        println(Crayon(foreground=:red, bold=true)("Performing RMsynthesis"))
        FDF, realFDF, imagFDF = RMSynthesis(Qnu, Unu, nuArray * 1e6, PhiArray; log_progress = log_progress)
        Pmax = maxCube(FDF)
        WriteData3D(resultspath, FDF, "FDF", PhiArray)
        WriteData3D(resultspath, realFDF, "realFDF", PhiArray)
        WriteData3D(resultspath, imagFDF, "imagFDF", PhiArray)
        WriteData2D(resultspath, Pmax, "Pmax")
    else
        println("No Faraday tomography performed.")
    end
end

function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String,
                       df::DataFrame, add_noise, Noise_nu, kernel_size_synchrotron, IonizationFraction::Float64,
                       nuArray::AbstractArray, PhiArray, PixelLength_pc, PixelLength_cm,
                       BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false)
    println("-------------------------------------------")
    println("Processing Synchrotron data for LOS: $LOS")
    resultspath = joinpath(simu, LOS, "Synchrotron")
    mkpath(resultspath)  
    
    # Read simulation parameters
    SimuParameters = ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    B1, B2, BLOS = SimuParameters[1], SimuParameters[2], SimuParameters[3]
    T, n = SimuParameters[4], SimuParameters[5]
    SimuParameters = nothing

    los_pixel_length_pc, los_pixel_length_cm, los_distance_array = compute_los_spacing(BoxLength_pc, size(B1, 3))
    PixelLength_pc, PixelLength_cm, DistanceArray =
        los_pixel_length_pc, los_pixel_length_cm, los_distance_array

    Bperpcube = Bperp(B1, B2)
    psi_src = IntrinsicAngle(B1, B2)

    faraday_enabled = uppercase(FaradayRotation) == "Y"
    filtering_enabled = uppercase(responseSynchrotron) == "Y"
    noise_enabled = uppercase(add_noise) == "Y"

    # Compute electron density using alternative prescription
    println("-------------------------------------------")
    println("Computing electron density")
    ne = ne_propto_nH(n, IonizationFraction)
    WriteData3D(resultspath, ne, "ne", los_distance_array)

    # Compute integral of quantities
    println("-------------------------------------------")
    println("Computing integrated quantities")
    Btotal = Btot(B1,B2,BLOS)
    intBtotal = intLOS(Btotal, los_pixel_length_cm)
    sigmaBtotal = sigmaLOS(Btotal)
    Btotal = nothing
    Ne = intLOS(ne, los_pixel_length_cm)
    sigmane = sigmaLOS(ne)
    sigmaT = sigmaLOS(T)
    intBLOS = intLOS(BLOS, los_pixel_length_cm)
    sigmaBLOS = sigmaLOS(BLOS)
    intBperp = intLOS(Bperpcube, los_pixel_length_cm)
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
    if faraday_enabled
        resultspath = joinpath(resultspath, "WithFaraday")
        mkpath(resultspath)
        println("-------------------------------------------")
        println("Computing RM")
        dRM = deltaRM(BLOS, ne, los_pixel_length_pc)
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
    if faraday_enabled
        println("-------------------------------------------")
        println("Computing Qnu and Unu with Faraday rotation")
        Qnu, Unu = QUnu3D(Bperpcube, psi_src, RMcube, nuArray, df, los_pixel_length_cm)
    else
        println("-------------------------------------------")
        println("Computing Qnu and Unu without Faraday rotation")
        Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, los_pixel_length_cm)
    end
    T_nu = Tnu3D(Bperpcube, nuArray, df, los_pixel_length_cm)
    Bperpcube = nothing
    
    if filtering_enabled
        kernel = Kernel.gaussian(kernel_size_synchrotron)
        println("-------------------------------------------")
        println("Applying filtering")
        @views for i in axes(Qnu, 3)
            Qnu[:, :, i] = HighPass(Qnu[:, :, i], kernel)
            Unu[:, :, i] = HighPass(Unu[:, :, i], kernel)
            T_nu[:, :, i] = HighPass(T_nu[:, :, i], kernel)
        end
        resultspath = joinpath(resultspath, "filtered")
        mkpath(resultspath)
    else
        println("No filtering performed.")
    end

    # Adding noise
    if noise_enabled
        println("-------------------------------------------")
        println("Adding a gaussian noise to Q and U with sigma = ", Noise_nu)
        noiseQ = similar(Qnu[:, :, 1])
        noiseU = similar(Unu[:, :, 1])
        distQ = Normal(0, Noise_nu)
        distU = Normal(0, Noise_nu)
        @views for i in axes(Qnu, 3)
            rand!(distQ, noiseQ)
            rand!(distU, noiseU)
            Qnu[:, :, i] .+= noiseQ
            Unu[:, :, i] .+= noiseU
        end
    end

    WriteData3D(resultspath, Qnu, "Qnu", nuArray)
    WriteData3D(resultspath, Unu, "Unu", nuArray)
    WriteData3D(resultspath, T_nu, "T_nu", nuArray)
    mv(joinpath(resultspath, "T_nu.fits"), joinpath(resultspath, "Tnu.fits"),force=true)  # Rename file
    
    println("-------------------------------------------")
    println("Computing Pnu and Pnumax")
    Pnucube = Pnu(Qnu, Unu)
    Pnumax = maxCube(Pnucube)
    WriteData3D(resultspath, Pnucube, "Pnu", nuArray)
    WriteData2D(resultspath, Pnumax, "Pnumax")
end


function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String,
                       df::DataFrame,  add_noise, Noise_nu, kernel_size_synchrotron, nuArray::AbstractArray, PhiArray,
                       PixelLength_pc, PixelLength_cm, BoxLength_pc,
                       DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false)

    println("-------------------------------------------")
    println("Processing Synchrotron data for LOS: $LOS")
    resultspath = joinpath(simu, LOS, "Synchrotron")
    mkpath(resultspath)  # Create directory if it doesn't exist
    
    # Read simulation parameters
    SimuParameters = ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    B1, B2, BLOS = SimuParameters[1], SimuParameters[2], SimuParameters[3]
    T, n = SimuParameters[4], SimuParameters[5]
    nHp = SimuParameters[7]
    SimuParameters = nothing

    los_pixel_length_pc, los_pixel_length_cm, los_distance_array = compute_los_spacing(BoxLength_pc, size(B1, 3))
    PixelLength_pc, PixelLength_cm, DistanceArray =
        los_pixel_length_pc, los_pixel_length_cm, los_distance_array

    Bperpcube = Bperp(B1, B2)
    psi_src = IntrinsicAngle(B1, B2)

    faraday_enabled = uppercase(FaradayRotation) == "Y"
    filtering_enabled = uppercase(responseSynchrotron) == "Y"
    noise_enabled = uppercase(add_noise) == "Y"

    # Compute integral of quantities
    println("-------------------------------------------")
    println("Computing integrated quantities")
    Btotal = Btot(B1,B2,BLOS)
    intBtotal = intLOS(Btotal, los_pixel_length_cm)
    sigmaBtotal = sigmaLOS(Btotal)
    Btotal = nothing
    Ne = intLOS(nHp, los_pixel_length_cm)
    sigmane = sigmaLOS(nHp)
    sigmaT = sigmaLOS(T)
    intBLOS = intLOS(BLOS, los_pixel_length_cm)
    sigmaBLOS = sigmaLOS(BLOS)
    intBperp = intLOS(Bperpcube, los_pixel_length_cm)
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
    if faraday_enabled
        resultspath = joinpath(resultspath, "WithFaraday")
        mkpath(resultspath)
        println("Computing RM")
        dRM = deltaRM(BLOS, nHp, los_pixel_length_pc)
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
    if faraday_enabled
        println("-------------------------------------------")
        println("Computing Qnu and Unu with Faraday rotation")
        Qnu, Unu = QUnu3D(Bperpcube, psi_src, RMcube, nuArray, df, los_pixel_length_cm)
    else
        println("-------------------------------------------")
        println("Computing Qnu and Unu without Faraday rotation")
        Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, los_pixel_length_cm)
    end
    T_nu = Tnu3D(Bperpcube, nuArray, df, los_pixel_length_cm)
    Bperpcube = nothing

    if filtering_enabled
        kernel = Kernel.gaussian(kernel_size_synchrotron)
        println("-------------------------------------------")
        println("Applying filtering")
        @views for i in axes(Qnu, 3)
            Qnu[:, :, i] = HighPass(Qnu[:, :, i], kernel)
            Unu[:, :, i] = HighPass(Unu[:, :, i], kernel)
            T_nu[:, :, i] = HighPass(T_nu[:, :, i], kernel)
        end
        resultspath = joinpath(resultspath, "filtered")
        mkpath(resultspath)
    else
        println("No filtering performed.")
    end
    # Adding noise
    if noise_enabled
        println("-------------------------------------------")
        println("Adding a gaussian noise to Q and U with sigma = ", Noise_nu)
        noiseQ = similar(Qnu[:, :, 1])
        noiseU = similar(Unu[:, :, 1])
        distQ = Normal(0, Noise_nu)
        distU = Normal(0, Noise_nu)
        @views for i in axes(Qnu, 3)
            rand!(distQ, noiseQ)
            rand!(distU, noiseU)
            Qnu[:, :, i] .+= noiseQ
            Unu[:, :, i] .+= noiseU
        end
    end

    WriteData3D(resultspath, Qnu, "Qnu", nuArray)
    WriteData3D(resultspath, Unu, "Unu", nuArray)
    WriteData3D(resultspath, T_nu, "T_nu", nuArray)
    mv(joinpath(resultspath, "T_nu.fits"), joinpath(resultspath, "Tnu.fits"),force=true)  # Rename file

    println("-------------------------------------------")
    println("Computing Pnu and Pnumax")
    Pnucube = Pnu(Qnu, Unu)
    Pnumax = maxCube(Pnucube)
    WriteData3D(resultspath, Pnucube, "Pnu", nuArray)
    WriteData2D(resultspath, Pnumax, "Pnumax")

    if faraday_enabled
        println("-------------------------------------------")
        println("Performing RMsynthesis")
        FDF, realFDF, imagFDF = RMSynthesis(Qnu, Unu, nuArray * 1e6, PhiArray; log_progress = log_progress)
        Pmax = maxCube(FDF)
        WriteData3D(resultspath, FDF, "FDF", PhiArray)
        WriteData3D(resultspath, realFDF, "realFDF", PhiArray)
        WriteData3D(resultspath, imagFDF, "imagFDF", PhiArray)
        WriteData2D(resultspath, Pmax, "Pmax")
    else
        println("No Faraday tomography performed.")
    end
end
