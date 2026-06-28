using DataFrames

function compute_los_spacing(box_length_pc, depth)
    depth > 0 || error("Cannot compute line-of-sight spacing for an empty cube.")
    los_pixel_length_pc = Float64(box_length_pc) / depth
    los_pixel_length_cm = los_pixel_length_pc * PARSEC_TO_CM
    los_distance_array = range(start=0.0, step=los_pixel_length_pc, length=depth)
    return los_pixel_length_pc, los_pixel_length_cm, los_distance_array
end

_stage(message) = @info message

function _validate_processing_cube_shapes(simu, LOS, expected_shape, cubes::Pair...)
    expected_shape === nothing && return nothing
    expected = Tuple(Int.(expected_shape))

    for (name, cube) in cubes
        cube === nothing && continue
        actual = size(cube)
        actual == expected || throw_config_error(
            "Cube shape mismatch for $(name) in $(simu) after LOS=$(LOS): expected $(expected) from `BoxLength_pix`, got $(actual). " *
            "Update `BoxLength_pix` or check the FITS cube dimensions.";
            code=:cube_shape_mismatch,
        )
    end

    return nothing
end

function _write_integrated_quantities(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm; metadata=nothing)
    Btotal = Btot(B1, B2, BLOS)
    intBtotal = intLOS(Btotal, los_pixel_length_cm)
    sigmaBtotal = sigmaLOS(Btotal)
    Ne = intLOS(ne, los_pixel_length_cm)
    sigmane = sigmaLOS(ne)
    sigmaT = sigmaLOS(T)
    intBLOS = intLOS(BLOS, los_pixel_length_cm)
    sigmaBLOS = sigmaLOS(BLOS)
    intBperp = intLOS(Bperpcube, los_pixel_length_cm)

    WriteData2D(resultspath, intBtotal, "intBtotal"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, sigmaBtotal, "sigmaBtotal"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, Ne, "intne"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, sigmane, "sigmane"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, sigmaT, "sigmaT"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, intBLOS, "intBLOS"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, sigmaBLOS, "sigmaBLOS"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, intBperp, "intBperp"; ensure_path=false, metadata=metadata)

    return nothing
end

function _healpix_unit(DataName::String)
    return haskey(DictHeader, DataName) ? String(DictHeader[DataName]["bunit"]) : ""
end

function _healpix_extname(DataName::String)
    return uppercase(replace(DataName, "_" => ""))
end

function _healpix_map_vector(data)
    if data isa AbstractVector
        return collect(data)
    elseif ndims(data) == 2 && size(data, 2) == 1
        return collect(view(data, :, 1))
    end

    error("Expected a HEALPix map vector or Npix x 1 array, got size $(size(data)).")
end

function _healpix_stack_matrix(data)
    ndims(data) == 3 && size(data, 2) == 1 ||
        error("Expected a HEALPix stack encoded as Npix x 1 x Nslice, got size $(size(data)).")
    return reshape(data, size(data, 1), size(data, 3))
end

function _write_healpix_map_quantity(resultspath, data, DataName::String, hp_meta)
    path = joinpath(resultspath, "$(DataName).fits")
    write_healpix_map(path, _healpix_map_vector(data);
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit(DataName),
        extname = _healpix_extname(DataName),
    )
    return path
end

function _write_healpix_stack_quantity(resultspath, data, DataName::String, coordinates, hp_meta; basename::Union{Nothing,String}=nothing)
    matrix = _healpix_stack_matrix(data)
    if size(matrix, 2) == 1
        return [_write_healpix_map_quantity(resultspath, view(matrix, :, 1), DataName, hp_meta)]
    end

    return write_healpix_stack(resultspath, matrix, basename === nothing ? DataName : basename, collect(coordinates);
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit(DataName),
        extname = _healpix_extname(DataName),
    )
end

function _write_healpix_rmclean_result(resultspath, result, hp_meta)
    write_healpix_stack(resultspath, Matrix(result.cleanFDF), "cleanFDF", result.phi;
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit("cleanFDF"),
        extname = _healpix_extname("cleanFDF"),
    )
    write_healpix_stack(resultspath, Matrix(result.realCleanFDF), "realCleanFDF", result.phi;
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit("realCleanFDF"),
        extname = _healpix_extname("realCleanFDF"),
    )
    write_healpix_stack(resultspath, Matrix(result.imagCleanFDF), "imagCleanFDF", result.phi;
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit("imagCleanFDF"),
        extname = _healpix_extname("imagCleanFDF"),
    )
    write_healpix_stack(resultspath, Matrix(abs.(result.residual)), "residualFDF", result.phi;
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit("residualFDF"),
        extname = _healpix_extname("residualFDF"),
    )

    return nothing
end

function _write_rmclean_cartesian(resultspath, result; metadata=nothing)
    WriteData3D(resultspath, result.cleanFDF, "cleanFDF", result.phi; ensure_path=false, metadata=metadata)
    WriteData3D(resultspath, result.realCleanFDF, "realCleanFDF", result.phi; ensure_path=false, metadata=metadata)
    WriteData3D(resultspath, result.imagCleanFDF, "imagCleanFDF", result.phi; ensure_path=false, metadata=metadata)
    WriteData3D(resultspath, abs.(result.residual), "residualFDF", result.phi; ensure_path=false, metadata=metadata)
    return nothing
end

function _write_integrated_quantities_healpix(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm, hp_meta)
    Btotal = Btot(B1, B2, BLOS)
    quantities = (
        "intBtotal" => intLOS(Btotal, los_pixel_length_cm),
        "sigmaBtotal" => sigmaLOS(Btotal),
        "intne" => intLOS(ne, los_pixel_length_cm),
        "sigmane" => sigmaLOS(ne),
        "sigmaT" => sigmaLOS(T),
        "intBLOS" => intLOS(BLOS, los_pixel_length_cm),
        "sigmaBLOS" => sigmaLOS(BLOS),
        "intBperp" => intLOS(Bperpcube, los_pixel_length_cm),
    )

    for (name, data) in quantities
        _write_healpix_map_quantity(resultspath, data, name, hp_meta)
    end

    return nothing
end

function _apply_synchrotron_filter!(Qnu, Unu, T_nu, Llarge_filter_pix)
    # The filter works entirely in PIXEL units (Δx = Δy = 1 pixel), matching
    # the documented convention: `kernel_size_synchrotron` is the largest
    # retained spatial scale in pixels. No small-scale cut is applied beyond
    # the Nyquist limit (Lcut_small = 2 pixels ⇔ f = fNy = 0.5 cycle/pixel).
    n, m = size(Qnu, 1), size(Qnu, 2)
    H, _ = instrument_bandpass_L(
        n,
        m;
        Δx = 1.0,
        Δy = 1.0,
        Lcut_small = 2.0,
        Llarge = Float64(Llarge_filter_pix),
        fNy = 0.5,
    )

    Qnu .= apply_to_array_xy(Qnu, H; n = n, m = m)
    Unu .= apply_to_array_xy(Unu, H; n = n, m = m)
    T_nu .= apply_to_array_xy(T_nu, H; n = n, m = m)

    return nothing
end

"""
    _add_noise!(Qnu, Unu, SNR_nu, rng)

Add gaussian noise to the Q and U cubes so that the per-channel polarized
signal-to-noise ratio equals `SNR_nu`. The reference signal in each frequency
channel is the rms of the polarized intensity, P_rms = sqrt(<Q²> + <U²>), and
the same standard deviation σ = P_rms / SNR_nu is applied to Q and U
(σ in the same unit as the cubes, i.e. Kelvin).
"""
function _add_noise!(Qnu, Unu, SNR_nu, rng)
    SNR_nu > 0 || error("SNR_nu must be > 0, got $SNR_nu")
    noiseQ = similar(Qnu[:, :, 1])
    noiseU = similar(Unu[:, :, 1])

    @views for i in axes(Qnu, 3)
        Qch = Qnu[:, :, i]
        Uch = Unu[:, :, i]
        P_rms = sqrt(mean(abs2, Qch) + mean(abs2, Uch))
        sigma = P_rms / SNR_nu
        sigma > 0 || continue
        randn!(rng, noiseQ)
        randn!(rng, noiseU)
        Qch .+= sigma .* noiseQ
        Uch .+= sigma .* noiseU
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
    SNR_nu,
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
    expected_shape = nothing,
    metadata = nothing,
    rm_clean_enabled::Bool = false,
    rm_clean_gain::Real = 0.1,
    rm_clean_niter::Integer = 1000,
    rm_clean_threshold::Real = 0.0,
)
    resultspath = joinpath(simu, LOS, "Synchrotron")
    mkpath(resultspath)
    fits_metadata = metadata === nothing ? Dict{String, Any}() : copy(metadata)
    fits_metadata["LOS"] = String(LOS)
    fits_metadata["BLENPC"] = Float64(BoxLength_pc)
    grid_kind = simulation_grid_kind(simu)
    hp_meta = grid_kind == :healpix ? healpix_simulation_metadata(simu) : nothing
    fits_metadata["GRID"] = grid_kind == :healpix ? "HEALPIX" : "IMAGE"

    SimuParameters = ReadSimulation(simu, LOS, conversionn, conversionT, conversionB)
    B1, B2, BLOS = SimuParameters[1], SimuParameters[2], SimuParameters[3]
    T, n = SimuParameters[4], SimuParameters[5]
    nHp = SimuParameters[7]
    _validate_processing_cube_shapes(simu, LOS, expected_shape,
        "B1" => B1,
        "B2" => B2,
        "BLOS" => BLOS,
        "temperature" => T,
        "density" => n,
        "densityHp" => nHp,
    )

    los_pixel_length_pc, los_pixel_length_cm, los_distance_array = compute_los_spacing(BoxLength_pc, size(B1, 3))

    Bperpcube = Bperp(B1, B2)
    psi_src = IntrinsicAngle(B1, B2)

    faraday_enabled = uppercase(FaradayRotation) == "Y"
    filtering_enabled = uppercase(responseSynchrotron) == "Y"
    noise_enabled = uppercase(add_noise) == "Y"

    _stage("Computing electron density")
    ne = electron_density_builder(T, n, nHp)
    if write_ne
        if hp_meta === nothing
            WriteData3D(resultspath, ne, "ne", los_distance_array; ensure_path=false, metadata=fits_metadata)
        else
            _write_healpix_stack_quantity(resultspath, ne, "ne", los_distance_array, hp_meta)
        end
    end

    _stage("Computing integrated quantities")
    if hp_meta === nothing
        _write_integrated_quantities(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm; metadata=fits_metadata)
    else
        _write_integrated_quantities_healpix(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm, hp_meta)
    end
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
        if hp_meta === nothing
            WriteData2D(resultspath, RMmap, "RMmap"; ensure_path=false, metadata=fits_metadata)
        else
            _write_healpix_map_quantity(resultspath, RMmap, "RMmap", hp_meta)
        end
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
        hp_meta === nothing ||
            throw_config_error("HEALPix inputs do not support the current cartesian Fourier filtering (`responseSynchrotron=Y`). Disable filtering or add a spherical-harmonic HEALPix filter.";
                code=:unsupported_grid_operation)
        _stage("Applying interferometric Fourier mask")
        _apply_synchrotron_filter!(Qnu, Unu, T_nu, kernel_size_synchrotron)
        resultspath = joinpath(resultspath, "filtered")
        mkpath(resultspath)
    else
        @info "No filtering performed"
    end

    if noise_enabled
        _stage("Adding gaussian noise to Q and U")
        _add_noise!(Qnu, Unu, SNR_nu, rng)
    end

    # The internal frequency axis is in MHz; FITS headers declare CUNIT3 = "Hz",
    # so the spectral axis is converted to Hz once, here, at write time.
    nuArray_Hz = collect(Float64, nuArray) .* 1e6

    if hp_meta === nothing
        WriteQUnu3D(resultspath, Qnu, Unu, nuArray_Hz; ensure_path=false, metadata=fits_metadata)
        WriteData3D(resultspath, T_nu, "T_nu", nuArray_Hz; ensure_path=false, metadata=fits_metadata, filename="Tnu.fits")
    else
        _write_healpix_stack_quantity(resultspath, Qnu, "Qnu", nuArray_Hz, hp_meta)
        _write_healpix_stack_quantity(resultspath, Unu, "Unu", nuArray_Hz, hp_meta)
        _write_healpix_stack_quantity(resultspath, T_nu, "T_nu", nuArray_Hz, hp_meta; basename="Tnu")
    end

    _stage("Computing Pnu and Pnumax")
    Pnucube = Pnu(Qnu, Unu)
    Pnumax = maxCube(Pnucube)
    if hp_meta === nothing
        WriteData3D(resultspath, Pnucube, "Pnu", nuArray_Hz; ensure_path=false, metadata=fits_metadata)
        WriteData2D(resultspath, Pnumax, "Pnumax"; ensure_path=false, metadata=fits_metadata)
    else
        _write_healpix_stack_quantity(resultspath, Pnucube, "Pnu", nuArray_Hz, hp_meta)
        _write_healpix_map_quantity(resultspath, Pnumax, "Pnumax", hp_meta)
    end
    try
        write_polarization_diagnostic_plots(resultspath, Qnu, Unu, T_nu, nuArray_Hz; Pnumax = Pnumax)
    catch err
        @warn "Failed to write polarization diagnostic plots" path = resultspath exception = (err, catch_backtrace())
    end

    if faraday_enabled
        _stage("Performing RM synthesis")
        rmsf = nothing
        if hp_meta === nothing
            FDF, realFDF, imagFDF = RMSynthesis(Qnu, Unu, nuArray * 1e6, PhiArray; log_progress = log_progress)
            Pmax = maxCube(FDF)
            WriteData3D(resultspath, FDF, "FDF", PhiArray; ensure_path=false, metadata=fits_metadata)
            WriteData3D(resultspath, realFDF, "realFDF", PhiArray; ensure_path=false, metadata=fits_metadata)
            WriteData3D(resultspath, imagFDF, "imagFDF", PhiArray; ensure_path=false, metadata=fits_metadata)
            WriteData2D(resultspath, Pmax, "Pmax"; ensure_path=false, metadata=fits_metadata)
        else
            q_stack = HealpixStack(_healpix_stack_matrix(Qnu); nside=hp_meta.nside, order=hp_meta.order)
            u_stack = HealpixStack(_healpix_stack_matrix(Unu); nside=hp_meta.nside, order=hp_meta.order)
            hp_result = RMSynthesisHealpix(q_stack, u_stack, nuArray * 1e6, PhiArray; log_progress = log_progress)
            write_healpix_rm_result(resultspath, hp_result; overwrite=true)
            Pmax = maximum(hp_result.fdf, dims=2)
            _write_healpix_map_quantity(resultspath, Pmax, "Pmax", hp_meta)
        end

        # RMSF diagnostics: write the spread function and its resolution metrics
        # alongside the FDF products. If RM-CLEAN is requested, diagnostics are
        # required and failures should abort the run; otherwise they only warn.
        if rm_clean_enabled
            rmsf = rmsf_diagnostics(nuArray * 1e6, PhiArray)
            @info "RMSF diagnostics" fwhm = rmsf.fwhm delta_phi_theory = rmsf.fwhm_theoretical phi_max = rmsf.phi_max max_scale = rmsf.max_scale
            write_rmsf(resultspath, rmsf; ensure_path=false, metadata=fits_metadata)
        else
            try
                rmsf = rmsf_diagnostics(nuArray * 1e6, PhiArray)
                @info "RMSF diagnostics" fwhm = rmsf.fwhm delta_phi_theory = rmsf.fwhm_theoretical phi_max = rmsf.phi_max max_scale = rmsf.max_scale
                write_rmsf(resultspath, rmsf; ensure_path=false, metadata=fits_metadata)
            catch err
                @warn "Failed to compute or write RMSF diagnostics" exception = err
            end
        end

        if rm_clean_enabled
            _stage("Performing RM-CLEAN")
            if hp_meta === nothing
                clean_result = rmclean(realFDF, imagFDF, PhiArray, rmsf;
                    gain = rm_clean_gain,
                    threshold = rm_clean_threshold,
                    niter = rm_clean_niter,
                    log_progress = log_progress,
                )
                _write_rmclean_cartesian(resultspath, clean_result; metadata=fits_metadata)
            else
                clean_result = RMClean(q_stack.pixels, u_stack.pixels, nuArray * 1e6, PhiArray;
                    gain = rm_clean_gain,
                    threshold = rm_clean_threshold,
                    niter = rm_clean_niter,
                    diagnostics = rmsf,
                    log_progress = log_progress,
                )
                _write_healpix_rmclean_result(resultspath, clean_result, hp_meta)
            end
        end
    else
        @info "No Faraday tomography performed"
    end

    return nothing
end

function ProcessSynchrotron(simu::AbstractString, LOS, FaradayRotation::AbstractString, responseSynchrotron::AbstractString,
                       df::DataFrame, add_noise, SNR_nu, kernel_size_synchrotron, zeta::Float64, Geff::Float64,
                       omegaPAH::Float64, XC::Float64, nuArray::AbstractArray, PhiArray,
                       PixelLength_pc, PixelLength_cm, BoxLength_pc,
                       DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng(), expected_shape = nothing, metadata = nothing,
                       rm_clean_enabled::Bool = false, rm_clean_gain::Real = 0.1, rm_clean_niter::Integer = 1000, rm_clean_threshold::Real = 0.0)
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        SNR_nu,
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
        expected_shape = expected_shape,
        metadata = metadata,
        rm_clean_enabled = rm_clean_enabled,
        rm_clean_gain = rm_clean_gain,
        rm_clean_niter = rm_clean_niter,
        rm_clean_threshold = rm_clean_threshold,
    )
end

function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String,
                       df::DataFrame, add_noise, SNR_nu, kernel_size_synchrotron, IonizationFraction::Float64,
                       nuArray::AbstractArray, PhiArray, PixelLength_pc, PixelLength_cm,
                       BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng(), expected_shape = nothing, metadata = nothing,
                       rm_clean_enabled::Bool = false, rm_clean_gain::Real = 0.1, rm_clean_niter::Integer = 1000, rm_clean_threshold::Real = 0.0)
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        SNR_nu,
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
        expected_shape = expected_shape,
        metadata = metadata,
        rm_clean_enabled = rm_clean_enabled,
        rm_clean_gain = rm_clean_gain,
        rm_clean_niter = rm_clean_niter,
        rm_clean_threshold = rm_clean_threshold,
    )
end

function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String,
                       df::DataFrame,  add_noise, SNR_nu, kernel_size_synchrotron, nuArray::AbstractArray, PhiArray,
                       PixelLength_pc, PixelLength_cm, BoxLength_pc,
                       DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng(), expected_shape = nothing, metadata = nothing,
                       rm_clean_enabled::Bool = false, rm_clean_gain::Real = 0.1, rm_clean_niter::Integer = 1000, rm_clean_threshold::Real = 0.0)
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        SNR_nu,
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
        expected_shape = expected_shape,
        metadata = metadata,
        rm_clean_enabled = rm_clean_enabled,
        rm_clean_gain = rm_clean_gain,
        rm_clean_niter = rm_clean_niter,
        rm_clean_threshold = rm_clean_threshold,
    )
end
