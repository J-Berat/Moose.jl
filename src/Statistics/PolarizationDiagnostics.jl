function _unwrap_phase(phases::AbstractVector)
    isempty(phases) && return Float64[]

    unwrapped = collect(Float64, phases)
    offset = 0.0
    previous = NaN
    @inbounds for i in eachindex(unwrapped)
        current = unwrapped[i]
        if !isfinite(current)
            previous = NaN
            offset = 0.0
            continue
        end
        if !isfinite(previous)
            previous = current
            continue
        end
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

const _latexstrings_pkgid = Base.PkgId(Base.UUID("b964fa9f-0449-5b57-a5c2-d3ea65f4040f"), "LaTeXStrings")

function _ltx(text::AbstractString)
    LaTeXStrings = Base.require(_latexstrings_pkgid)
    return Base.invokelatest(getproperty(LaTeXStrings, :LaTeXString), text)
end

function _finite_extrema(values::AbstractVector)
    finite_values = [Float64(value) for value in values if isfinite(value)]
    isempty(finite_values) && return (-1.0, 1.0)
    return extrema(finite_values)
end

function _padded_limits(values::AbstractVector; pad_fraction::Real = 0.08)
    lo, hi = _finite_extrema(values)
    span = hi - lo
    if !isfinite(span) || span <= 0
        width = max(abs(lo), 1.0) * 0.1
        return (lo - width, hi + width)
    end

    pad = span * Float64(pad_fraction)
    return (lo - pad, hi + pad)
end

function _equal_span_limits(xvalues::AbstractVector, yvalues::AbstractVector; pad_fraction::Real = 0.10)
    xlo, xhi = _padded_limits(xvalues; pad_fraction = pad_fraction)
    ylo, yhi = _padded_limits(yvalues; pad_fraction = pad_fraction)
    xmid = (xlo + xhi) / 2
    ymid = (ylo + yhi) / 2
    halfspan = max(xhi - xlo, yhi - ylo) / 2
    halfspan = isfinite(halfspan) && halfspan > 0 ? halfspan : 1.0
    return (xmid - halfspan, xmid + halfspan), (ymid - halfspan, ymid + halfspan)
end

function _zero_line_if_visible!(CairoMakie, ax, xlims, ylims)
    ylims[1] <= 0 <= ylims[2] &&
        _cm(CairoMakie, Symbol("hlines!"), ax, [0.0]; color = (:gray45, 0.35), linewidth = 1)
    xlims[1] <= 0 <= xlims[2] &&
        _cm(CairoMakie, Symbol("vlines!"), ax, [0.0]; color = (:gray45, 0.35), linewidth = 1)
    return nothing
end

function _representative_stokes_spectra(Qnu::AbstractArray, Unu::AbstractArray, Tnu::AbstractArray, Pnumax)
    size(Qnu) == size(Unu) == size(Tnu) || error("Qnu, Unu and Tnu must have the same shape.")
    ndims(Qnu) >= 1 || error("Qnu, Unu and Tnu must include at least one spectral dimension.")

    if ndims(Qnu) == 1
        return collect(Float64, Qnu), collect(Float64, Unu), collect(Float64, Tnu), ()
    end

    if Pnumax === nothing
        spatial_index = Tuple(fill(1, ndims(Qnu) - 1))
    else
        finite_index = findfirst(isfinite, Pnumax)
        finite_index === nothing && error("Pnumax contains no finite pixel from which to select a representative spectrum.")
        spatial_index = Tuple(finite_index)
        best = Pnumax[finite_index]
        for idx in eachindex(Pnumax)
            value = Pnumax[idx]
            if isfinite(value) && value > best
                best = value
                spatial_index = Tuple(CartesianIndices(Pnumax)[idx])
            end
        end
    end
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
    all(nu -> isfinite(nu) && nu > 0, nu_hz) ||
        throw(ArgumentError("Frequencies must be positive and finite."))

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

function _draw_polarization_angle!(CairoMakie, ax, spectra)
    _cm(CairoMakie, Symbol("scatterlines!"), ax, spectra.lambda2, spectra.psi_deg;
        color = :black, markersize = 6, linewidth = 2)
    xlims = _padded_limits(spectra.lambda2)
    ylims = _padded_limits(spectra.psi_deg)
    _zero_line_if_visible!(CairoMakie, ax, xlims, ylims)
    _cm(CairoMakie, Symbol("limits!"), ax, xlims[1], xlims[2], ylims[1], ylims[2])
    return nothing
end

function _draw_fractional_polarization!(CairoMakie, ax, spectra; legend::Bool = true)
    _cm(CairoMakie, Symbol("scatterlines!"), ax, spectra.lambda2, spectra.frac_p;
        color = :black, markersize = 6, linewidth = 2, label = _ltx(raw"\mathrm{Total}\;p"))
    _cm(CairoMakie, Symbol("scatterlines!"), ax, spectra.lambda2, spectra.frac_q;
        color = :blue, markersize = 5, linewidth = 2, label = _ltx(raw"\mathrm{Stokes}\;q"))
    _cm(CairoMakie, Symbol("scatterlines!"), ax, spectra.lambda2, spectra.frac_u;
        color = :red, markersize = 5, linewidth = 2, label = _ltx(raw"\mathrm{Stokes}\;u"))
    xlims = _padded_limits(spectra.lambda2)
    ylims = _padded_limits(vcat(spectra.frac_p, spectra.frac_q, spectra.frac_u))
    _zero_line_if_visible!(CairoMakie, ax, xlims, ylims)
    legend && _cm(CairoMakie, :axislegend, ax; position = :rb, labelsize = 20)
    _cm(CairoMakie, Symbol("limits!"), ax, xlims[1], xlims[2], ylims[1], ylims[2])
    return nothing
end

function _draw_stokes_qu!(CairoMakie, fig_cell, colorbar_cell, spectra; colorbar::Bool = true)
    ax = _cm(CairoMakie, :Axis, fig_cell;
        xlabel = _ltx(raw"\mathrm{Stokes}\;q"),
        ylabel = _ltx(raw"\mathrm{Stokes}\;u"),
        title = _ltx(raw"\mathrm{Stokes}\;q\!-\!u\;\mathrm{trajectory}"),
        aspect = _cm(CairoMakie, :DataAspect))
    _cm(CairoMakie, Symbol("lines!"), ax, spectra.frac_q, spectra.frac_u;
        color = (:black, 0.55), linewidth = 1.5)
    scatter = _cm(CairoMakie, Symbol("scatter!"), ax, spectra.frac_q, spectra.frac_u;
        color = spectra.lambda2,
        colormap = :coolwarm,
        strokecolor = :black,
        strokewidth = 0.45,
        markersize = 7)
    xlims, ylims = _equal_span_limits(spectra.frac_q, spectra.frac_u)
    _zero_line_if_visible!(CairoMakie, ax, xlims, ylims)
    if colorbar
        _cm(CairoMakie, :Colorbar, colorbar_cell, scatter;
            label = _ltx(raw"\lambda^2\;(\mathrm{m}^2)"),
            labelsize = 22,
            ticklabelsize = 20,
            width = 18)
    end
    _cm(CairoMakie, Symbol("limits!"), ax, xlims[1], xlims[2], ylims[1], ylims[2])
    return ax
end

function _save_polarization_angle_plot(path::AbstractString, spectra)
    CairoMakie = _ensure_cairomakie()
    fig = _cm(CairoMakie, :Figure; size = (980, 620), fontsize = 24)
    ax = _cm(CairoMakie, :Axis, _figure_cell(fig, 1, 1);
        xlabel = _ltx(raw"\lambda^2\;(\mathrm{m}^2)"),
        ylabel = _ltx(raw"\Psi\;(\mathrm{degrees})"),
        title = _ltx(raw"\mathrm{Polarization\ angle}"))
    _draw_polarization_angle!(CairoMakie, ax, spectra)
    atomic_write_path(path) do tmp_path
        _cm(CairoMakie, :save, tmp_path, fig)
    end
end

function _save_fractional_polarization_plot(path::AbstractString, spectra)
    CairoMakie = _ensure_cairomakie()
    fig = _cm(CairoMakie, :Figure; size = (980, 620), fontsize = 24)
    ax = _cm(CairoMakie, :Axis, _figure_cell(fig, 1, 1);
        xlabel = _ltx(raw"\lambda^2\;(\mathrm{m}^2)"),
        ylabel = _ltx(raw"\mathrm{Fractional\ polarization}"),
        title = _ltx(raw"\mathrm{Fractional\ Stokes\ parameters}"))
    _draw_fractional_polarization!(CairoMakie, ax, spectra)
    atomic_write_path(path) do tmp_path
        _cm(CairoMakie, :save, tmp_path, fig)
    end
end

function _save_stokes_qu_plot(path::AbstractString, spectra)
    CairoMakie = _ensure_cairomakie()
    fig = _cm(CairoMakie, :Figure; size = (820, 760), fontsize = 24)
    _draw_stokes_qu!(CairoMakie, _figure_cell(fig, 1, 1), _figure_cell(fig, 1, 2), spectra)
    atomic_write_path(path) do tmp_path
        _cm(CairoMakie, :save, tmp_path, fig)
    end
end

function _save_polarization_composite_plot(path::AbstractString, spectra)
    CairoMakie = _ensure_cairomakie()
    fig = _cm(CairoMakie, :Figure; size = (1900, 620), fontsize = 22)

    ax_angle = _cm(CairoMakie, :Axis, _figure_cell(fig, 1, 1);
        xlabel = _ltx(raw"\lambda^2\;(\mathrm{m}^2)"),
        ylabel = _ltx(raw"\Psi\;(\mathrm{degrees})"),
        title = _ltx(raw"\mathrm{Polarization\ angle}"))
    _draw_polarization_angle!(CairoMakie, ax_angle, spectra)

    ax_fraction = _cm(CairoMakie, :Axis, _figure_cell(fig, 1, 2);
        xlabel = _ltx(raw"\lambda^2\;(\mathrm{m}^2)"),
        ylabel = _ltx(raw"\mathrm{Fractional\ polarization}"),
        title = _ltx(raw"\mathrm{Fractional\ Stokes\ parameters}"))
    _draw_fractional_polarization!(CairoMakie, ax_fraction, spectra; legend = true)

    _draw_stokes_qu!(CairoMakie, _figure_cell(fig, 1, 3), _figure_cell(fig, 1, 4), spectra)

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
        composite_png = joinpath(resultspath, "polarization_diagnostics.png"),
        composite_pdf = joinpath(resultspath, "polarization_diagnostics.pdf"),
    )

    _save_polarization_angle_plot(paths.angle, spectra)
    _save_fractional_polarization_plot(paths.fraction, spectra)
    _save_stokes_qu_plot(paths.qu, spectra)
    _save_polarization_composite_plot(paths.composite_png, spectra)
    _save_polarization_composite_plot(paths.composite_pdf, spectra)

    @info "Wrote polarization diagnostic plots" path = resultspath pixel = spectra.spatial_index
    return paths
end
