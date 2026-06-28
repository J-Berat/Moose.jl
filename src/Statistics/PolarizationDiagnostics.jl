function _unwrap_phase(phases::AbstractVector)
    isempty(phases) && return Float64[]

    unwrapped = collect(Float64, phases)
    offset = 0.0
    previous = unwrapped[1]
    @inbounds for i in 2:length(unwrapped)
        current = unwrapped[i]
        delta = current - previous
        if delta > pi
            offset -= 2pi
        elseif delta < -pi
            offset += 2pi
        end
        previous = current
        unwrapped[i] = current + offset
    end

    return unwrapped
end

function _safe_fraction(numerator::Real, denominator::Real)
    isfinite(numerator) && isfinite(denominator) && abs(denominator) > eps(Float64) || return NaN
    return Float64(numerator) / Float64(denominator)
end

function _cm(CairoMakie, name::Symbol, args...; kwargs...)
    return Base.invokelatest(getproperty(CairoMakie, name), args...; kwargs...)
end

_figure_cell(fig, args...) = Base.invokelatest(getindex, fig, args...)

function _representative_stokes_spectra(Qnu::AbstractArray, Unu::AbstractArray, Tnu::AbstractArray, Pnumax)
    size(Qnu) == size(Unu) == size(Tnu) || error("Qnu, Unu and Tnu must have the same shape.")
    ndims(Qnu) >= 1 || error("Qnu, Unu and Tnu must include at least one spectral dimension.")

    if ndims(Qnu) == 1
        return collect(Float64, Qnu), collect(Float64, Unu), collect(Float64, Tnu), ()
    end

    spatial_index = Pnumax === nothing ? Tuple(fill(1, ndims(Qnu) - 1)) : Tuple(argmax(Pnumax))
    length(spatial_index) == ndims(Qnu) - 1 ||
        error("Pnumax shape $(size(Pnumax)) is incompatible with Stokes cube shape $(size(Qnu)).")

    selectors = (spatial_index..., :)
    return (
        collect(Float64, @view Qnu[selectors...]),
        collect(Float64, @view Unu[selectors...]),
        collect(Float64, @view Tnu[selectors...]),
        spatial_index,
    )
end

function polarization_diagnostic_spectra(Qnu::AbstractArray, Unu::AbstractArray, Tnu::AbstractArray,
    nu_hz::AbstractArray; Pnumax = nothing)

    q_stokes, u_stokes, i_stokes, spatial_index = _representative_stokes_spectra(Qnu, Unu, Tnu, Pnumax)
    length(q_stokes) == length(nu_hz) ||
        error("Frequency axis length $(length(nu_hz)) does not match Stokes spectral length $(length(q_stokes)).")

    lambda2 = (C_m ./ collect(Float64, nu_hz)) .^ 2
    p_stokes = hypot.(q_stokes, u_stokes)
    frac_q = [_safe_fraction(q_stokes[i], i_stokes[i]) for i in eachindex(q_stokes)]
    frac_u = [_safe_fraction(u_stokes[i], i_stokes[i]) for i in eachindex(u_stokes)]
    frac_p = [_safe_fraction(p_stokes[i], i_stokes[i]) for i in eachindex(p_stokes)]
    psi_deg = rad2deg.(0.5 .* _unwrap_phase(atan.(u_stokes, q_stokes)))

    return (;
        lambda2,
        psi_deg,
        frac_p,
        frac_q,
        frac_u,
        spatial_index,
    )
end

function _save_polarization_angle_plot(path::AbstractString, spectra)
    CairoMakie = _ensure_cairomakie()
    fig = _cm(CairoMakie, :Figure; size = (900, 620))
    ax = _cm(CairoMakie, :Axis, _figure_cell(fig, 1, 1);
        xlabel = "lambda^2 (m^2)",
        ylabel = "Psi (degrees)",
        title = "Polarization angle")
    _cm(CairoMakie, Symbol("scatterlines!"), ax, spectra.lambda2, spectra.psi_deg; color = :black, markersize = 7)
    _cm(CairoMakie, Symbol("hlines!"), ax, [0.0]; color = (:gray45, 0.45), linewidth = 1)
    _cm(CairoMakie, Symbol("tightlimits!"), ax)
    atomic_write_path(path) do tmp_path
        _cm(CairoMakie, :save, tmp_path, fig)
    end
end

function _save_fractional_polarization_plot(path::AbstractString, spectra)
    CairoMakie = _ensure_cairomakie()
    fig = _cm(CairoMakie, :Figure; size = (900, 620))
    ax = _cm(CairoMakie, :Axis, _figure_cell(fig, 1, 1);
        xlabel = "lambda^2 (m^2)",
        ylabel = "Fractional polarization",
        title = "Fractional Stokes parameters")
    _cm(CairoMakie, Symbol("scatterlines!"), ax, spectra.lambda2, spectra.frac_p; color = :black, markersize = 7, label = "Total p")
    _cm(CairoMakie, Symbol("scatterlines!"), ax, spectra.lambda2, spectra.frac_q; color = :blue, markersize = 6, label = "Stokes q")
    _cm(CairoMakie, Symbol("scatterlines!"), ax, spectra.lambda2, spectra.frac_u; color = :red, markersize = 6, label = "Stokes u")
    _cm(CairoMakie, Symbol("hlines!"), ax, [0.0]; color = (:gray45, 0.45), linewidth = 1)
    _cm(CairoMakie, :axislegend, ax; position = :rb)
    _cm(CairoMakie, Symbol("tightlimits!"), ax)
    atomic_write_path(path) do tmp_path
        _cm(CairoMakie, :save, tmp_path, fig)
    end
end

function _save_stokes_qu_plot(path::AbstractString, spectra)
    CairoMakie = _ensure_cairomakie()
    fig = _cm(CairoMakie, :Figure; size = (760, 620))
    ax = _cm(CairoMakie, :Axis, _figure_cell(fig, 1, 1);
        xlabel = "Stokes q",
        ylabel = "Stokes u",
        title = "Stokes q-u trajectory",
        aspect = _cm(CairoMakie, :DataAspect))
    scatter = _cm(CairoMakie, Symbol("scatter!"), ax, spectra.frac_q, spectra.frac_u;
        color = spectra.lambda2,
        colormap = :coolwarm,
        strokecolor = :black,
        strokewidth = 1,
        markersize = 12)
    _cm(CairoMakie, Symbol("hlines!"), ax, [0.0]; color = (:gray45, 0.5), linewidth = 1)
    _cm(CairoMakie, Symbol("vlines!"), ax, [0.0]; color = (:gray45, 0.5), linewidth = 1)
    _cm(CairoMakie, :Colorbar, _figure_cell(fig, 1, 2), scatter; label = "lambda^2 (m^2)")
    _cm(CairoMakie, Symbol("tightlimits!"), ax)
    atomic_write_path(path) do tmp_path
        _cm(CairoMakie, :save, tmp_path, fig)
    end
end

function write_polarization_diagnostic_plots(resultspath::AbstractString, Qnu::AbstractArray, Unu::AbstractArray,
    Tnu::AbstractArray, nu_hz::AbstractArray; Pnumax = nothing)

    spectra = polarization_diagnostic_spectra(Qnu, Unu, Tnu, nu_hz; Pnumax = Pnumax)
    paths = (;
        angle = joinpath(resultspath, "polarization_angle_vs_lambda2.png"),
        fraction = joinpath(resultspath, "fractional_polarization_vs_lambda2.png"),
        qu = joinpath(resultspath, "stokes_qu_diagram.png"),
    )

    _save_polarization_angle_plot(paths.angle, spectra)
    _save_fractional_polarization_plot(paths.fraction, spectra)
    _save_stokes_qu_plot(paths.qu, spectra)

    @info "Wrote polarization diagnostic plots" path = resultspath pixel = spectra.spatial_index
    return paths
end
