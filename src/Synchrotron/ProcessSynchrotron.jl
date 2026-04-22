using DataFrames

function compute_los_spacing(box_length_pc, depth)
    depth > 0 || error("Cannot compute line-of-sight spacing for an empty cube.")
    los_pixel_length_pc = Float64(box_length_pc) / depth
    los_pixel_length_cm = los_pixel_length_pc * PARSEC_TO_CM
    los_distance_array = range(start=0.0, step=los_pixel_length_pc, length=depth)
    return los_pixel_length_pc, los_pixel_length_cm, los_distance_array
end

_stage(message) = @info message

function _write_integrated_quantities(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm)
    Btotal = Btot(B1, B2, BLOS)
    intBtotal = intLOS(Btotal, los_pixel_length_cm)
    sigmaBtotal = sigmaLOS(Btotal)
    Ne = intLOS(ne, los_pixel_length_cm)
    sigmane = sigmaLOS(ne)
    sigmaT = sigmaLOS(T)
    intBLOS = intLOS(BLOS, los_pixel_length_cm)
    sigmaBLOS = sigmaLOS(BLOS)
    intBperp = intLOS(Bperpcube, los_pixel_length_cm)

    WriteData2D(resultspath, intBtotal, "intBtotal"; ensure_path=false)
    WriteData2D(resultspath, sigmaBtotal, "sigmaBtotal"; ensure_path=false)
    WriteData2D(resultspath, Ne, "intne"; ensure_path=false)
    WriteData2D(resultspath, sigmane, "sigmane"; ensure_path=false)
    WriteData2D(resultspath, sigmaT, "sigmaT"; ensure_path=false)
    WriteData2D(resultspath, intBLOS, "intBLOS"; ensure_path=false)
    WriteData2D(resultspath, sigmaBLOS, "sigmaBLOS"; ensure_path=false)
    WriteData2D(resultspath, intBperp, "intBperp"; ensure_path=false)

    return nothing
end

function _apply_synchrotron_filter!(Qnu, Unu, T_nu, Llarge_filter, box_length_pc)
    n, m = size(Qnu, 1), size(Qnu, 2)
    Δx = Float64(box_length_pc) / n
    Δy = Float64(box_length_pc) / m
    fNy = min(1 / (2Δx), 1 / (2Δy))
    H, _ = instrument_bandpass_L(
        n,
        m;
        Δx = Δx,
        Δy = Δy,
        Lcut_small = 1.0,
        Llarge = Float64(Llarge_filter),
        fNy = fNy,
    )

    Qnu .= apply_to_array_xy(Qnu, H; n = n, m = m)
    Unu .= apply_to_array_xy(Unu, H; n = n, m = m)
    T_nu .= apply_to_array_xy(T_nu, H; n = n, m = m)

    return nothing
end

function _add_noise!(Qnu, Unu, Noise_nu, rng)
    noiseQ = similar(Qnu[:, :, 1])
    noiseU = similar(Unu[:, :, 1])
    distQ = Normal(0, Noise_nu)
    distU = Normal(0, Noise_nu)

    @views for i in axes(Qnu, 3)
        rand!(rng, distQ, noiseQ)
        rand!(rng, distU, noiseU)
        Qnu[:, :, i] .+= noiseQ
        Unu[:, :, i] .+= noiseU
    end

    return nothing
end

function _process_synchrotron_common(
    simu::AbstractString,
    LOS,
    FaradayRotation::AbstractString,
    responseSynchrotron::AbstractString,
    df::DataFrame,
    add_noise,
    Noise_nu,
    kernel_size_synchrotron,
    nuArray::AbstractArray,
    PhiArray,
    BoxLength_pc,
    conversionn,
    conversionT,
    conversionB,
    electron_density_builder;
    write_ne::Bool = true,
    log_progress::Bool = false,
    rng = Random.default_rng(),
)
    resultspath = joinpath(simu, LOS, "Synchrotron")
    mkpath(resultspath)

    SimuParameters = ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    B1, B2, BLOS = SimuParameters[1], SimuParameters[2], SimuParameters[3]
    T, n = SimuParameters[4], SimuParameters[5]
    nHp = SimuParameters[7]

    los_pixel_length_pc, los_pixel_length_cm, los_distance_array = compute_los_spacing(BoxLength_pc, size(B1, 3))

    Bperpcube = Bperp(B1, B2)
    psi_src = IntrinsicAngle(B1, B2)

    faraday_enabled = uppercase(FaradayRotation) == "Y"
    filtering_enabled = uppercase(responseSynchrotron) == "Y"
    noise_enabled = uppercase(add_noise) == "Y"

    _stage("Computing electron density")
    ne = electron_density_builder(T, n, nHp)
    write_ne && WriteData3D(resultspath, ne, "ne", los_distance_array; ensure_path=false)

    _stage("Computing integrated quantities")
    _write_integrated_quantities(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm)
    B1 = nothing
    B2 = nothing
    T = nothing
    n = nothing

    RMcube = nothing
    if faraday_enabled
        resultspath = joinpath(resultspath, "WithFaraday")
        mkpath(resultspath)
        _stage("Computing RM")
        dRM = deltaRM(BLOS, ne, los_pixel_length_pc)
        RMcube = RM(dRM)
        RMmap = RMcube[:, :, end]
        WriteData2D(resultspath, RMmap, "RMmap"; ensure_path=false)
    else
        resultspath = joinpath(resultspath, "noFaraday")
        mkpath(resultspath)
        @info "No Faraday rotation included"
    end

    if faraday_enabled
        _stage("Computing Qnu and Unu with Faraday rotation")
        Qnu, Unu = QUnu3D(Bperpcube, psi_src, RMcube, nuArray, df, los_pixel_length_cm; log_progress = log_progress)
    else
        _stage("Computing Qnu and Unu without Faraday rotation")
        Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, los_pixel_length_cm; log_progress = log_progress)
    end

    T_nu = Tnu3D(Bperpcube, nuArray, df, los_pixel_length_cm)
    Bperpcube = nothing

    if filtering_enabled
        _stage("Applying interferometric Fourier mask")
        _apply_synchrotron_filter!(Qnu, Unu, T_nu, kernel_size_synchrotron, BoxLength_pc)
        resultspath = joinpath(resultspath, "filtered")
        mkpath(resultspath)
    else
        @info "No filtering performed"
    end

    if noise_enabled
        _stage("Adding gaussian noise to Q and U")
        _add_noise!(Qnu, Unu, Noise_nu, rng)
    end

    WriteQUnu3D(resultspath, Qnu, Unu, nuArray; ensure_path=false)
    WriteData3D(resultspath, T_nu, "T_nu", nuArray; ensure_path=false)
    mv(joinpath(resultspath, "T_nu.fits"), joinpath(resultspath, "Tnu.fits"), force=true)

    _stage("Computing Pnu and Pnumax")
    Pnucube = Pnu(Qnu, Unu)
    Pnumax = maxCube(Pnucube)
    WriteData3D(resultspath, Pnucube, "Pnu", nuArray; ensure_path=false)
    WriteData2D(resultspath, Pnumax, "Pnumax"; ensure_path=false)

    if faraday_enabled
        _stage("Performing RM synthesis")
        FDF, realFDF, imagFDF = RMSynthesis(Qnu, Unu, nuArray * 1e6, PhiArray; log_progress = log_progress)
        Pmax = maxCube(FDF)
        WriteData3D(resultspath, FDF, "FDF", PhiArray; ensure_path=false)
        WriteData3D(resultspath, realFDF, "realFDF", PhiArray; ensure_path=false)
        WriteData3D(resultspath, imagFDF, "imagFDF", PhiArray; ensure_path=false)
        WriteData2D(resultspath, Pmax, "Pmax"; ensure_path=false)
    else
        @info "No Faraday tomography performed"
    end

    return nothing
end

function ProcessSynchrotron(simu::AbstractString, LOS, FaradayRotation::AbstractString, responseSynchrotron::AbstractString,
                       df::DataFrame, add_noise, Noise_nu, kernel_size_synchrotron, zeta::Float64, Geff::Float64,
                       omegaPAH::Float64, XC::Float64, nuArray::AbstractArray, PhiArray,
                       PixelLength_pc, PixelLength_cm, BoxLength_pc,
                       DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng())
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        Noise_nu,
        kernel_size_synchrotron,
        nuArray,
        PhiArray,
        BoxLength_pc,
        conversionn,
        conversionT,
        conversionB,
        (T, n, nHp) -> Wolfire_ne(zeta, Geff, omegaPAH, XC, T, n);
        write_ne = true,
        log_progress = log_progress,
        rng = rng,
    )
end

function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String,
                       df::DataFrame, add_noise, Noise_nu, kernel_size_synchrotron, IonizationFraction::Float64,
                       nuArray::AbstractArray, PhiArray, PixelLength_pc, PixelLength_cm,
                       BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng())
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        Noise_nu,
        kernel_size_synchrotron,
        nuArray,
        PhiArray,
        BoxLength_pc,
        conversionn,
        conversionT,
        conversionB,
        (T, n, nHp) -> ne_propto_nH(n, IonizationFraction);
        write_ne = true,
        log_progress = log_progress,
        rng = rng,
    )
end

function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String,
                       df::DataFrame,  add_noise, Noise_nu, kernel_size_synchrotron, nuArray::AbstractArray, PhiArray,
                       PixelLength_pc, PixelLength_cm, BoxLength_pc,
                       DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng())
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        Noise_nu,
        kernel_size_synchrotron,
        nuArray,
        PhiArray,
        BoxLength_pc,
        conversionn,
        conversionT,
        conversionB,
        (T, n, nHp) -> nHp;
        write_ne = false,
        log_progress = log_progress,
        rng = rng,
    )
end
