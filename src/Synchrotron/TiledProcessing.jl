# Tiled (banded) synchrotron processing for inputs larger than RAM.
#
# When `tile_size` is set, the sky is processed in bands. For cartesian cubes
# the bands are `tile_size` rows along the second sky axis; for HEALPix inputs
# they are `tile_size` HEALPix pixel rows of the Npix x Nshell table. Input
# data are read from FITS one band at a time, all per-pixel products are
# computed on the band, 3D outputs are streamed into pre-created FITS files,
# and 2D maps are assembled in memory (a sky map is small compared to the
# cubes). Every per-pixel computation is identical to the non-tiled path, so
# tiled and non-tiled runs produce the same values.
#
# Not supported in tiled mode (guarded in the config validation and here):
# interferometric filtering (needs the full sky plane in Fourier space),
# noise injection (the per-channel σ is derived from the full-map rms),
# RM-CLEAN, and the polarization diagnostic plots (which need the full Q/U
# cubes; they are skipped with an @info message). Tiled HEALPix runs stream
# their 3D products as single-file HEALPix cubes (see `write_healpix_cube`)
# instead of one FITS file per slice, and require each simulation field to be
# a single-column HEALPix binary table.

"""
    _band_ranges(n, band) -> Vector{UnitRange{Int}}

Split `1:n` into consecutive ranges of length `band` (last one possibly
shorter).
"""
function _band_ranges(n::Integer, band::Integer)
    band > 0 || error("Band size must be positive, got $band.")
    return [i:min(i + band - 1, n) for i in 1:band:n]
end

# Permutation applied by `permute_dims(A, LOS)` (FITSUtils): original cube
# axes -> processing axes (sky1, sky2, LOS).
_los_perm(LOS::AbstractString) = LOS == "x" ? (2, 3, 1) : LOS == "y" ? (3, 1, 2) : (1, 2, 3)

# Shape of the processing-frame cube given the on-disk shape.
_permuted_shape(sz::NTuple{3, Int}, LOS) = ntuple(d -> sz[_los_perm(LOS)[d]], 3)

"""
    _read_band(hdu, LOS, jrange, conversion, ::Type{T}) -> Array{T, 3}

Read the sky band `[:, jrange, :]` of the processing-frame cube from an
on-disk cube, applying the same LOS axis permutation as `permute_dims`, the
unit conversion, and the working precision. The returned band has shape
`(n1, length(jrange), n_los)`.
"""
function _read_band(hdu, LOS::AbstractString, jrange::UnitRange{Int}, conversion, ::Type{T}) where {T <: AbstractFloat}
    raw = if LOS == "z"
        read(hdu, :, jrange, :)
    elseif LOS == "x"
        # processing frame (a2, a3, a1): band rows live on original axis 3.
        permutedims(read(hdu, :, :, jrange), (2, 3, 1))
    elseif LOS == "y"
        # processing frame (a3, a1, a2): band rows live on original axis 1.
        permutedims(read(hdu, jrange, :, :), (3, 1, 2))
    else
        error("Unknown LOS: $LOS (expected \"x\", \"y\" or \"z\")")
    end

    band = T.(raw .* conversion)
    bad = findfirst(x -> !isfinite(x), band)
    bad === nothing || error("FITS cube band contains a non-finite value at index $(bad) (band rows $(jrange)).")
    return band
end

function _open_tiled_cube(file)
    fits = FITS(file, "r")
    for i in 1:length(fits)
        hdu = fits[i]
        if hdu isa ImageHDU && ndims(hdu) == 3
            return fits, hdu
        end
    end
    close(fits)
    error("No 3D image HDU found in $(file); `tile_size` requires plain 3D FITS cubes.")
end

function _tiled_field_handle(simu, field; field_sources=nothing)
    source = simulation_field_source(simu, field, field_sources)
    source isa AbstractString ||
        throw_config_error("`tile_size` requires each simulation field to be a single FITS cube; $(field) resolves to multiple files.";
            code=:invalid_tile_size)
    validate_required_fits(source)
    return _open_tiled_cube(source)
end

"""
    _StreamedCube3D

FITS cube written band by band. The file is created up front with its final
dimensions and full header, filled with direct `ffppxll` writes (one
contiguous block per channel and band), and atomically renamed on close.
"""
mutable struct _StreamedCube3D
    fits::FITS
    final_path::String
    tmp_path::String
    dims::NTuple{3, Int}
end

function _streamed_cube(resultspath, DataName::String, dims::NTuple{3, Int}, specarray, ::Type{T};
    metadata = nothing, filename = nothing) where {T <: AbstractFloat}
    params = _header_params_cached(DataName)
    header = buildHeader3D(
        params["naxis"], dims,
        params["ctype1"], params["ctype2"], params["ctype3"],
        params["cunit1"], params["cunit2"], params["cunit3"],
        params["bunit"], specarray;
        metadata = metadata,
    )

    fits_name = filename === nothing ? "$DataName.fits" : String(filename)
    final_path = joinpath(resultspath, fits_name)
    tmp_path = final_path * ".part"
    rm(tmp_path; force = true)

    f = FITS(tmp_path, "w")
    FITSIO.fits_create_img(f.fitsfile, T, Int64[dims...])
    FITSIO.fits_write_header(f.fitsfile, header, true)
    return _StreamedCube3D(f, final_path, tmp_path, dims)
end

# cfitsio shared library, via FITSIO's own import of CFITSIO (both are pinned
# in the Manifest, so no new dependency is introduced).
const _LIBCFITSIO = FITSIO.libcfitsio

"""
    _fits_write_packed(fitsfile, fpixel, nelements, buffer_ptr, ::Type{T})

Write `nelements` consecutive on-disk pixels starting at the FITS **file**
coordinates `fpixel`, reading the values from the packed buffer behind
`buffer_ptr` (element type `T`). This calls cfitsio's `ffppxll`
(`fits_write_pix`) directly: the high-level Julia wrappers bounds-check
`data` as if it were the full-size image, so they cannot stream a small band
block at an offset file position. The caller must `GC.@preserve` the array
that backs `buffer_ptr`.
"""
function _fits_write_packed(fitsfile, fpixel::Vector{Int64}, nelements::Integer, buffer_ptr::Ptr, ::Type{T}) where {T}
    status = Ref{Cint}(0)
    ccall(
        (:ffppxll, _LIBCFITSIO),
        Cint,
        (Ptr{Cvoid}, Cint, Ptr{Int64}, Int64, Ptr{Cvoid}, Ref{Cint}),
        fitsfile.ptr,
        FITSIO.CFITSIO.cfitsio_typecode(T),
        fpixel,
        Int64(nelements),
        buffer_ptr,
        status,
    )
    FITSIO.fits_assert_ok(status[])
    return nothing
end

function _write_band!(sc::_StreamedCube3D, block::AbstractArray{<:Real, 3}, jrange::UnitRange{Int})
    n1, _, nch = sc.dims
    size(block) == (n1, length(jrange), nch) ||
        error("Band block has shape $(size(block)); expected ($(n1), $(length(jrange)), $(nch)).")

    # For each channel, the band occupies one contiguous run of the on-disk
    # (column-major) pixel stream, starting at file pixel (1, first(jrange), k);
    # the matching source values are the k-th slice of the packed band block.
    data = block isa Array ? block : Array(block)
    nel = n1 * length(jrange)
    GC.@preserve data begin
        for k in 1:nch
            _fits_write_packed(
                sc.fits.fitsfile,
                Int64[1, first(jrange), k],
                nel,
                pointer(data, (k - 1) * nel + 1),
                eltype(data),
            )
        end
    end
    return nothing
end

function _finalize!(sc::_StreamedCube3D)
    close(sc.fits)
    mv(sc.tmp_path, sc.final_path; force = true)
    @info "Wrote FITS file" path = sc.final_path
    return nothing
end

function _abort!(sc::_StreamedCube3D)
    try
        close(sc.fits)
    catch
    end
    rm(sc.tmp_path; force = true)
    return nothing
end

# ---------------------------------------------------------------------------
# HEALPix tiled I/O
# ---------------------------------------------------------------------------

const _CF = FITSIO.CFITSIO

# Open handle on a single-column HEALPix binary table, positioned on the
# table HDU for partial-row reads.
struct _TiledHealpixHandle
    fits::FITS
    colnum::Int
    npix::Int
    nshell::Int
end

function _open_tiled_healpix(file::AbstractString)
    fits = FITS(String(file), "r")
    info = nothing
    try
        info = _healpix_table_hdu_info(fits)
        info === nothing && error("No HEALPix binary table HDU found in $(file).")
        if length(info.colnames) != 1
            throw_config_error(
                "`tile_size` with HEALPix inputs requires a single-column HEALPix table; $(file) has $(length(info.colnames)) columns ($(join(info.colnames, ", "))).";
                code=:invalid_tile_size)
        end
        npix = Healpix.nside2npix(info.nside)
        info.nrows == npix ||
            error("HEALPix table in $(file) has $(info.nrows) rows, expected $(npix) for NSIDE=$(info.nside).")

        _CF.fits_movabs_hdu(fits.fitsfile, info.hdu_index)
        _, repeat, _ = _CF.fits_get_coltype(fits.fitsfile, 1)
        return _TiledHealpixHandle(fits, 1, npix, Int(repeat)), info
    catch
        close(fits)
        rethrow()
    end
end

"""
    _read_healpix_band(handle, jr, conversion) -> Array{Float64, 3}

Read HEALPix pixel rows `jr` from a tiled input handle as an
`(nrows, 1, nshell)` processing cube, converting UNSEEN sentinels to `NaN`
and applying the unit conversion.
"""
function _read_healpix_band(handle::_TiledHealpixHandle, jr::UnitRange{Int}, conversion)
    nrows = length(jr)
    buffer = Vector{Float64}(undef, handle.nshell * nrows)
    _CF.fits_read_col(handle.fits.fitsfile, handle.colnum, first(jr), 1, buffer)

    band = permutedims(reshape(buffer, handle.nshell, nrows))  # nrows x nshell
    _mask_unseen!(band)
    conversion == 1 || (band .*= conversion)
    return reshape(band, nrows, 1, handle.nshell)
end

"""
    _StreamedHealpixCube

Single-file HEALPix cube (see [`write_healpix_cube`](@ref)) written band by
band: the binary table is created up front with its final number of rows,
filled with partial `fits_write_col` writes, completed with a `COORDS`
extension, and atomically renamed on close.
"""
mutable struct _StreamedHealpixCube
    fits::FITS
    final_path::String
    tmp_path::String
    npix::Int
    nslice::Int
    coordinates::Vector{Float64}
    coordname::String
end

function _streamed_healpix_cube(resultspath, DataName::String, npix::Int, coordinates, hp_meta;
    filename = nothing, coordname::AbstractString = "COORD")

    nslice = length(coordinates)
    fits_name = filename === nothing ? "$(DataName).fits" : String(filename)
    final_path = joinpath(resultspath, fits_name)
    tmp_path = final_path * ".part"
    rm(tmp_path; force = true)

    f = FITS(tmp_path, "w")
    _CF.fits_create_binary_tbl(f.fitsfile, npix,
        [("PIXELS", "$(nslice)D", _healpix_unit(DataName))], _healpix_extname(DataName))
    _CF.fits_update_key(f.fitsfile, "PIXTYPE", "HEALPIX", "HEALPix pixelisation")
    _CF.fits_update_key(f.fitsfile, "ORDERING", hp_meta.order == :nested ? "NESTED" : "RING", "Pixel ordering scheme")
    _CF.fits_update_key(f.fitsfile, "NSIDE", Int(hp_meta.nside), "HEALPix NSIDE")
    _CF.fits_update_key(f.fitsfile, "FIRSTPIX", 0, "First pixel index")
    _CF.fits_update_key(f.fitsfile, "LASTPIX", npix - 1, "Last pixel index")
    _CF.fits_update_key(f.fitsfile, "INDXSCHM", "IMPLICIT", "Pixel indexing scheme")
    coordsys = get(hp_meta, :coordsys, nothing)
    coordsys === nothing || _CF.fits_update_key(f.fitsfile, "COORDSYS", String(coordsys), "Pixelisation coordinate system")

    return _StreamedHealpixCube(f, final_path, tmp_path, npix, nslice, Float64.(collect(coordinates)), String(coordname))
end

function _write_healpix_band!(sc::_StreamedHealpixCube, block::AbstractArray{<:Real, 3}, jr::UnitRange{Int})
    size(block) == (length(jr), 1, sc.nslice) ||
        error("HEALPix band block has shape $(size(block)); expected ($(length(jr)), 1, $(sc.nslice)).")

    data = permutedims(reshape(Float64.(block), length(jr), sc.nslice))  # nslice x nrows
    @inbounds for i in eachindex(data)
        isnan(data[i]) && (data[i] = HEALPIX_UNSEEN)
    end
    _CF.fits_write_col(sc.fits.fitsfile, 1, first(jr), 1, vec(data))
    return nothing
end

function _finalize!(sc::_StreamedHealpixCube)
    _CF.fits_create_binary_tbl(sc.fits.fitsfile, sc.nslice, [(sc.coordname, "1D", "")], "COORDS")
    _CF.fits_update_key(sc.fits.fitsfile, "COORDNAM", sc.coordname, "Physical meaning of the slice coordinate")
    _CF.fits_write_col(sc.fits.fitsfile, 1, 1, 1, sc.coordinates)
    close(sc.fits)
    mv(sc.tmp_path, sc.final_path; force = true)
    @info "Wrote FITS file" path = sc.final_path
    return nothing
end

function _abort!(sc::_StreamedHealpixCube)
    try
        close(sc.fits)
    catch
    end
    rm(sc.tmp_path; force = true)
    return nothing
end

function _process_synchrotron_tiled_healpix(
    simu::AbstractString,
    LOS,
    FaradayRotation::AbstractString,
    df::DataFrame,
    nuArray::AbstractArray,
    PhiArray,
    BoxLength_pc,
    conversionn,
    conversionT,
    conversionB,
    electron_density_builder;
    write_ne::Bool,
    log_progress::Bool,
    tile_rows::Int,
    resultspath::AbstractString,
    hp_meta,
    field_sources=nothing,
    physical_mask=nothing,
    density_kind::AbstractString="number_density",
    mean_molecular_weight::Real=1.0,
    hydrogen_mass_g::Real=M_p,
)
    faraday_enabled = uppercase(FaradayRotation) == "Y"

    handles = Dict{String, _TiledHealpixHandle}()
    streams = _StreamedHealpixCube[]

    try
        fields = ["Bx", "By", "Bz", "density", "temperature"]
        nHp_file = simulation_field_source(simu, "densityHp", field_sources)
        has_nHp = nHp_file isa AbstractString && isfile(nHp_file)
        has_nHp && push!(fields, "densityHp")

        for field in fields
            source = simulation_field_source(simu, field, field_sources)
            source isa AbstractString || throw_config_error(
                "`tile_size` with HEALPix inputs requires each simulation field to be a single FITS file; $(field) resolves to multiple files.";
                code=:invalid_tile_size)
            validate_required_fits(source)
            handle, info = _open_tiled_healpix(source)
            if info.nside != hp_meta.nside || info.order != hp_meta.order
                close(handle.fits)
                throw_config_error(
                    "`tile_size` with HEALPix inputs requires all fields on the same grid; $(field) is NSIDE=$(info.nside)/$(info.order), expected NSIDE=$(hp_meta.nside)/$(hp_meta.order).";
                    code=:invalid_tile_size)
            end
            handles[field] = handle
        end

        npix = handles["Bx"].npix
        nshell = handles["Bx"].nshell
        for (field, handle) in handles
            (handle.npix == npix && handle.nshell == nshell) || throw_config_error(
                "HEALPix stack shape mismatch for $(field) in $(simu): expected $(npix) x $(nshell), got $(handle.npix) x $(handle.nshell).";
                code=:cube_shape_mismatch)
        end

        los_pixel_length_pc, los_pixel_length_cm, los_distance_array = compute_los_spacing(BoxLength_pc, nshell)
        nuArray_Hz = collect(Float64, nuArray) .* 1e6
        Nfreq = length(nuArray)

        if faraday_enabled
            resultspath = joinpath(resultspath, "WithFaraday")
        else
            resultspath = joinpath(resultspath, "noFaraday")
            @info "No Faraday rotation included"
        end
        mkpath(resultspath)
        root_resultspath = dirname(resultspath)

        # 2D maps assembled across bands (a sky map is cheap to hold).
        maps2d = Dict{String, Vector{Float64}}(
            name => Vector{Float64}(undef, npix) for name in
            ("intBtotal", "sigmaBtotal", "intne", "sigmane", "sigmaT",
             "intBLOS", "sigmaBLOS", "intBperp", "Pnumax", "polfracmax", "alpha", "alpha_err")
        )
        faraday_enabled && (maps2d["RMmap"] = Vector{Float64}(undef, npix))
        faraday_enabled && (maps2d["Pmax"] = Vector{Float64}(undef, npix))

        # Streamed 3D outputs (single-file HEALPix cubes).
        ne_stream = nothing
        if write_ne
            ne_stream = _streamed_healpix_cube(root_resultspath, "ne", npix, los_distance_array, hp_meta)
            push!(streams, ne_stream)
        end
        q_stream = _streamed_healpix_cube(resultspath, "Qnu", npix, nuArray_Hz, hp_meta)
        u_stream = _streamed_healpix_cube(resultspath, "Unu", npix, nuArray_Hz, hp_meta)
        t_stream = _streamed_healpix_cube(resultspath, "T_nu", npix, nuArray_Hz, hp_meta; filename = "Tnu.fits")
        p_stream = _streamed_healpix_cube(resultspath, "Pnu", npix, nuArray_Hz, hp_meta)
        polfrac_stream = _streamed_healpix_cube(resultspath, "polfrac", npix, nuArray_Hz, hp_meta)
        append!(streams, (q_stream, u_stream, t_stream, p_stream, polfrac_stream))

        fdf_streams = nothing
        if faraday_enabled
            fdf_streams = (
                _streamed_healpix_cube(resultspath, "FDF", npix, PhiArray, hp_meta; coordname = "PHI"),
                _streamed_healpix_cube(resultspath, "realFDF", npix, PhiArray, hp_meta; coordname = "PHI"),
                _streamed_healpix_cube(resultspath, "imagFDF", npix, PhiArray, hp_meta; coordname = "PHI"),
            )
            append!(streams, fdf_streams)
        end

        bands = _band_ranges(npix, tile_rows)
        for (band_index, jr) in enumerate(bands)
            _stage("Processing HEALPix pixel band $(band_index)/$(length(bands)) (pixels $(jr))")

            bx = _read_healpix_band(handles["Bx"], jr, conversionB)
            by = _read_healpix_band(handles["By"], jr, conversionB)
            bz = _read_healpix_band(handles["Bz"], jr, conversionB)
            B1, B2, BLOS = los_basis(bx, by, bz, LOS)
            T = _read_healpix_band(handles["temperature"], jr, conversionT)
            n = _read_healpix_band(handles["density"], jr, conversionn)
            nHp = has_nHp ? _read_healpix_band(handles["densityHp"], jr, conversionn) : nothing
            _density_to_number_density!(n, density_kind, mean_molecular_weight, hydrogen_mass_g)
            _apply_physical_mask!(B1, B2, BLOS, T, n, nHp, physical_mask)

            ne = electron_density_builder(T, n, nHp)
            ne_stream === nothing || _write_healpix_band!(ne_stream, ne, jr)

            Bperpband = Bperp(B1, B2)
            psi_src = IntrinsicAngle(B1, B2)
            Btotal = Btot(B1, B2, BLOS)

            maps2d["intBtotal"][jr] = vec(intLOS(Btotal, los_pixel_length_cm))
            maps2d["sigmaBtotal"][jr] = vec(sigmaLOS(Btotal))
            maps2d["intne"][jr] = vec(intLOS(ne, los_pixel_length_cm))
            maps2d["sigmane"][jr] = vec(sigmaLOS(ne))
            maps2d["sigmaT"][jr] = vec(sigmaLOS(T))
            maps2d["intBLOS"][jr] = vec(intLOS(BLOS, los_pixel_length_cm))
            maps2d["sigmaBLOS"][jr] = vec(sigmaLOS(BLOS))
            maps2d["intBperp"][jr] = vec(intLOS(Bperpband, los_pixel_length_cm))
            Btotal = nothing
            B1 = nothing
            B2 = nothing
            T = nothing
            n = nothing

            RMband = nothing
            if faraday_enabled
                dRM = deltaRM(BLOS, ne, los_pixel_length_pc)
                RMband = RM(dRM)
                maps2d["RMmap"][jr] = vec(RMband[:, :, end])
            end
            BLOS = nothing
            ne = nothing

            Qband, Uband = if faraday_enabled
                QUnu3D(Bperpband, psi_src, RMband, nuArray, df, los_pixel_length_cm; log_progress = false)
            else
                QUnuNoFaraday3D(Bperpband, psi_src, nuArray, df, los_pixel_length_cm; log_progress = false)
            end
            Tband = Tnu3D(Bperpband, nuArray, df, los_pixel_length_cm)
            Bperpband = nothing
            psi_src = nothing
            RMband = nothing

            _write_healpix_band!(q_stream, Qband, jr)
            _write_healpix_band!(u_stream, Uband, jr)
            _write_healpix_band!(t_stream, Tband, jr)

            Pband = Pnu(Qband, Uband)
            maps2d["Pnumax"][jr] = vec(maxCube(Pband))
            polfracband = PolarizationFraction(Pband, Tband)
            maps2d["polfracmax"][jr] = vec(_max_finite_cube(polfracband))
            _write_healpix_band!(p_stream, Pband, jr)
            _write_healpix_band!(polfrac_stream, polfracband, jr)
            Pband = nothing
            polfracband = nothing

            if Nfreq >= 2
                beta, alpha_err = spectral_index_map(Tband, nuArray; min_channels = 2)
                maps2d["alpha"][jr] = vec(beta) .+ 2.0
                maps2d["alpha_err"][jr] = vec(alpha_err)
            else
                maps2d["alpha"][jr] .= NaN
                maps2d["alpha_err"][jr] .= NaN
            end
            Tband = nothing

            if faraday_enabled
                FDF, realFDF, imagFDF = RMSynthesis(Qband, Uband, nuArray * 1e6, PhiArray; log_progress = false)
                maps2d["Pmax"][jr] = vec(maxCube(FDF))
                _write_healpix_band!(fdf_streams[1], FDF, jr)
                _write_healpix_band!(fdf_streams[2], realFDF, jr)
                _write_healpix_band!(fdf_streams[3], imagFDF, jr)
            end

            log_progress && print_progress(band_index, length(bands))
        end

        for sc in streams
            _finalize!(sc)
        end
        empty!(streams)

        # 2D products.
        for name in ("intBtotal", "sigmaBtotal", "intne", "sigmane", "sigmaT", "intBLOS", "sigmaBLOS", "intBperp")
            _write_healpix_map_quantity(root_resultspath, maps2d[name], name, hp_meta)
        end
        faraday_enabled && _write_healpix_map_quantity(resultspath, maps2d["RMmap"], "RMmap", hp_meta)
        _write_healpix_map_quantity(resultspath, maps2d["Pnumax"], "Pnumax", hp_meta)
        _write_healpix_map_quantity(resultspath, maps2d["polfracmax"], "polfracmax", hp_meta)
        faraday_enabled && _write_healpix_map_quantity(resultspath, maps2d["Pmax"], "Pmax", hp_meta)
        _write_healpix_map_quantity(resultspath, maps2d["alpha"], "alpha", hp_meta)
        _write_healpix_map_quantity(resultspath, maps2d["alpha_err"], "alpha_err", hp_meta)

        if faraday_enabled
            try
                rmsf = rmsf_diagnostics(nuArray * 1e6, PhiArray)
                @info "RMSF diagnostics" fwhm = rmsf.fwhm delta_phi_theory = rmsf.fwhm_theoretical phi_max = rmsf.phi_max max_scale = rmsf.max_scale
                write_rmsf(resultspath, rmsf; ensure_path = false)
            catch err
                @warn "Failed to compute or write RMSF diagnostics" exception = err
            end
        else
            @info "No Faraday tomography performed"
        end

        @info "Skipping polarization diagnostic plots in tiled mode (they require the full Q/U cubes in memory)"
    finally
        for sc in streams
            _abort!(sc)
        end
        for (_, handle) in handles
            try
                close(handle.fits)
            catch
            end
        end
    end

    return nothing
end

function _process_synchrotron_tiled(
    simu::AbstractString,
    LOS,
    FaradayRotation::AbstractString,
    df::DataFrame,
    nuArray::AbstractArray,
    PhiArray,
    BoxLength_pc,
    conversionn,
    conversionT,
    conversionB,
    electron_density_builder;
    write_ne::Bool,
    log_progress::Bool,
    expected_shape,
    fits_metadata,
    float_type::Type{<:AbstractFloat},
    tile_rows::Int,
    resultspath::AbstractString,
    field_sources=nothing,
    physical_mask=nothing,
    density_kind::AbstractString="number_density",
    mean_molecular_weight::Real=1.0,
    hydrogen_mass_g::Real=M_p,
)
    faraday_enabled = uppercase(FaradayRotation) == "Y"
    fits_metadata = copy(fits_metadata)
    fits_metadata["TILESIZE"] = tile_rows

    handles = Dict{String, Tuple{FITS, Any}}()
    streams = _StreamedCube3D[]

    try
        for field in ("Bx", "By", "Bz", "density", "temperature")
            handles[field] = _tiled_field_handle(simu, field; field_sources = field_sources)
        end
        nHp_file = simulation_field_source(simu, "densityHp", field_sources)
        has_nHp = nHp_file isa AbstractString && isfile(nHp_file)
        if has_nHp
            handles["densityHp"] = _open_tiled_cube(nHp_file)
        end

        disk_shape = NTuple{3, Int}(size(handles["Bx"][2]))
        for (field, (_, hdu)) in handles
            NTuple{3, Int}(size(hdu)) == disk_shape || throw_config_error(
                "Cube shape mismatch for $(field) in $(simu): expected $(disk_shape), got $(Tuple(size(hdu))).";
                code=:cube_shape_mismatch)
        end

        n1, n2, nz = _permuted_shape(disk_shape, LOS)
        if expected_shape !== nothing
            expected = Tuple(Int.(expected_shape))
            (n1, n2, nz) == expected || throw_config_error(
                "Cube shape mismatch in $(simu) after LOS=$(LOS): expected $(expected) from `BoxLength_pix`, got $((n1, n2, nz)). " *
                "Update `BoxLength_pix` or check the FITS cube dimensions.";
                code=:cube_shape_mismatch)
        end

        los_pixel_length_pc, los_pixel_length_cm, los_distance_array = compute_los_spacing(BoxLength_pc, nz)
        Nfreq = length(nuArray)
        nuArray_Hz = collect(Float64, nuArray) .* 1e6

        if faraday_enabled
            resultspath = joinpath(resultspath, "WithFaraday")
        else
            resultspath = joinpath(resultspath, "noFaraday")
            @info "No Faraday rotation included"
        end
        mkpath(resultspath)
        root_resultspath = dirname(resultspath)

        # 2D maps assembled across bands (a sky map is cheap to hold).
        maps2d = Dict{String, Matrix{Float64}}(
            name => Matrix{Float64}(undef, n1, n2) for name in
            ("intBtotal", "sigmaBtotal", "intne", "sigmane", "sigmaT",
             "intBLOS", "sigmaBLOS", "intBperp", "Pnumax", "polfracmax", "alpha", "alpha_err")
        )
        faraday_enabled && (maps2d["RMmap"] = Matrix{Float64}(undef, n1, n2))
        faraday_enabled && (maps2d["Pmax"] = Matrix{Float64}(undef, n1, n2))

        # Streamed 3D outputs.
        ne_stream = nothing
        if write_ne
            ne_stream = _streamed_cube(root_resultspath, "ne", (n1, n2, nz), los_distance_array, float_type; metadata = fits_metadata)
            push!(streams, ne_stream)
        end
        q_stream = _streamed_cube(resultspath, "Qnu", (n1, n2, Nfreq), nuArray_Hz, float_type; metadata = fits_metadata)
        u_stream = _streamed_cube(resultspath, "Unu", (n1, n2, Nfreq), nuArray_Hz, float_type; metadata = fits_metadata)
        t_stream = _streamed_cube(resultspath, "T_nu", (n1, n2, Nfreq), nuArray_Hz, float_type; metadata = fits_metadata, filename = "Tnu.fits")
        p_stream = _streamed_cube(resultspath, "Pnu", (n1, n2, Nfreq), nuArray_Hz, float_type; metadata = fits_metadata)
        polfrac_stream = _streamed_cube(resultspath, "polfrac", (n1, n2, Nfreq), nuArray_Hz, float_type; metadata = fits_metadata)
        append!(streams, (q_stream, u_stream, t_stream, p_stream, polfrac_stream))

        fdf_streams = nothing
        if faraday_enabled
            nPhi = length(PhiArray)
            fdf_streams = (
                _streamed_cube(resultspath, "FDF", (n1, n2, nPhi), PhiArray, float_type; metadata = fits_metadata),
                _streamed_cube(resultspath, "realFDF", (n1, n2, nPhi), PhiArray, float_type; metadata = fits_metadata),
                _streamed_cube(resultspath, "imagFDF", (n1, n2, nPhi), PhiArray, float_type; metadata = fits_metadata),
            )
            append!(streams, fdf_streams)
        end

        bands = _band_ranges(n2, tile_rows)
        for (band_index, jr) in enumerate(bands)
            _stage("Processing sky band $(band_index)/$(length(bands)) (rows $(jr))")

            bx = _read_band(handles["Bx"][2], LOS, jr, conversionB, float_type)
            by = _read_band(handles["By"][2], LOS, jr, conversionB, float_type)
            bz = _read_band(handles["Bz"][2], LOS, jr, conversionB, float_type)
            B1, B2, BLOS = los_basis(bx, by, bz, LOS)
            T = _read_band(handles["temperature"][2], LOS, jr, conversionT, float_type)
            n = _read_band(handles["density"][2], LOS, jr, conversionn, float_type)
            nHp = has_nHp ? _read_band(handles["densityHp"][2], LOS, jr, conversionn, float_type) : nothing
            _density_to_number_density!(n, density_kind, mean_molecular_weight, hydrogen_mass_g)
            _apply_physical_mask!(B1, B2, BLOS, T, n, nHp, physical_mask)

            ne = _to_precision(float_type, electron_density_builder(T, n, nHp))
            ne_stream === nothing || _write_band!(ne_stream, ne, jr)

            Bperpband = Bperp(B1, B2)
            psi_src = IntrinsicAngle(B1, B2)
            Btotal = Btot(B1, B2, BLOS)

            maps2d["intBtotal"][:, jr] = intLOS(Btotal, los_pixel_length_cm)
            maps2d["sigmaBtotal"][:, jr] = sigmaLOS(Btotal)
            maps2d["intne"][:, jr] = intLOS(ne, los_pixel_length_cm)
            maps2d["sigmane"][:, jr] = sigmaLOS(ne)
            maps2d["sigmaT"][:, jr] = sigmaLOS(T)
            maps2d["intBLOS"][:, jr] = intLOS(BLOS, los_pixel_length_cm)
            maps2d["sigmaBLOS"][:, jr] = sigmaLOS(BLOS)
            maps2d["intBperp"][:, jr] = intLOS(Bperpband, los_pixel_length_cm)
            Btotal = nothing
            B1 = nothing
            B2 = nothing
            T = nothing
            n = nothing

            RMband = nothing
            if faraday_enabled
                dRM = deltaRM(BLOS, ne, los_pixel_length_pc)
                RMband = RM(_to_precision(float_type, dRM))
                maps2d["RMmap"][:, jr] = RMband[:, :, end]
            end
            BLOS = nothing
            ne = nothing

            Qband, Uband = if faraday_enabled
                QUnu3D(Bperpband, psi_src, RMband, nuArray, df, los_pixel_length_cm; log_progress = false)
            else
                QUnuNoFaraday3D(Bperpband, psi_src, nuArray, df, los_pixel_length_cm; log_progress = false)
            end
            Tband = Tnu3D(Bperpband, nuArray, df, los_pixel_length_cm)
            Bperpband = nothing
            psi_src = nothing
            RMband = nothing

            _write_band!(q_stream, Qband, jr)
            _write_band!(u_stream, Uband, jr)
            _write_band!(t_stream, Tband, jr)

            Pband = Pnu(Qband, Uband)
            maps2d["Pnumax"][:, jr] = maxCube(Pband)
            polfracband = PolarizationFraction(Pband, Tband)
            maps2d["polfracmax"][:, jr] = _max_finite_cube(polfracband)
            _write_band!(p_stream, Pband, jr)
            _write_band!(polfrac_stream, polfracband, jr)
            Pband = nothing
            polfracband = nothing

            if Nfreq >= 2
                beta, alpha_err = spectral_index_map(Tband, nuArray; min_channels = 2)
                maps2d["alpha"][:, jr] = beta .+ 2.0
                maps2d["alpha_err"][:, jr] = alpha_err
            else
                maps2d["alpha"][:, jr] .= NaN
                maps2d["alpha_err"][:, jr] .= NaN
            end
            Tband = nothing

            if faraday_enabled
                FDF, realFDF, imagFDF = RMSynthesis(Qband, Uband, nuArray * 1e6, PhiArray; log_progress = false)
                maps2d["Pmax"][:, jr] = maxCube(FDF)
                _write_band!(fdf_streams[1], FDF, jr)
                _write_band!(fdf_streams[2], realFDF, jr)
                _write_band!(fdf_streams[3], imagFDF, jr)
            end

            log_progress && print_progress(band_index, length(bands))
        end

        for sc in streams
            _finalize!(sc)
        end
        empty!(streams)

        # 2D products.
        for name in ("intBtotal", "sigmaBtotal", "intne", "sigmane", "sigmaT", "intBLOS", "sigmaBLOS", "intBperp")
            WriteData2D(root_resultspath, maps2d[name], name; ensure_path = false, metadata = fits_metadata)
        end
        faraday_enabled && WriteData2D(resultspath, maps2d["RMmap"], "RMmap"; ensure_path = false, metadata = fits_metadata)
        WriteData2D(resultspath, maps2d["Pnumax"], "Pnumax"; ensure_path = false, metadata = fits_metadata)
        WriteData2D(resultspath, maps2d["polfracmax"], "polfracmax"; ensure_path = false, metadata = fits_metadata)
        faraday_enabled && WriteData2D(resultspath, maps2d["Pmax"], "Pmax"; ensure_path = false, metadata = fits_metadata)

        alpha_metadata = copy(fits_metadata)
        alpha_metadata["ALPHADEF"] = "S_nu ~ nu^alpha; alpha = beta_Tb + 2"
        WriteData2D(resultspath, maps2d["alpha"], "alpha"; ensure_path = false, metadata = alpha_metadata)
        WriteData2D(resultspath, maps2d["alpha_err"], "alpha_err"; ensure_path = false, metadata = alpha_metadata)

        if faraday_enabled
            try
                rmsf = rmsf_diagnostics(nuArray * 1e6, PhiArray)
                @info "RMSF diagnostics" fwhm = rmsf.fwhm delta_phi_theory = rmsf.fwhm_theoretical phi_max = rmsf.phi_max max_scale = rmsf.max_scale
                write_rmsf(resultspath, rmsf; ensure_path = false, metadata = fits_metadata)
            catch err
                @warn "Failed to compute or write RMSF diagnostics" exception = err
            end
        else
            @info "No Faraday tomography performed"
        end

        @info "Skipping polarization diagnostic plots in tiled mode (they require the full Q/U cubes in memory)"
    finally
        for sc in streams
            _abort!(sc)
        end
        for (_, (fits, _)) in handles
            try
                close(fits)
            catch
            end
        end
    end

    return nothing
end
