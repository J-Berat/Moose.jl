function plot_histogram(datasets, labels, x_label, y_label, bins, file_path, xlims=(nothing, nothing), ylims=(nothing, nothing))
    f = Figure(size = (1200, 800))
    ax = Axis(f[1,1], xlabel = x_label, ylabel = y_label, xgridvisible = false, ygridvisible = false)
    for (data, label) in zip(datasets, labels)
        hist!(ax, vec(filter(!isnan, data)), bins=bins, label=label)
        vlines!(ax, [mean(filter(!isnan, data))], label="Mean $label", linestyle=:dash)
    end
    xlims!(ax, xlims)
    ylims!(ax, ylims)
    f[1,2] = Legend(f, ax, framevisible=false)
    save(file_path, f)
end

function hist2D(
    x_data,
    y_data;
    T_values = nothing, 
    n_bins = 500,
    title = "2D Histogram with CairoMakie (log-log)",
    xlabel = "log10(n)",
    ylabel = "log10(P)",
    colorbar_label = "Counts (log scale)",
    percentage = 0.35,
)
    if any(x_data .<= 0)
        error("x_data contains non-positive values. Logarithmic transformation requires positive values.")
    end
    if any(y_data .<= 0)
        error("y_data contains non-positive values. Logarithmic transformation requires positive values.")
    end

    x_data = log10.(x_data)
    y_data = log10.(y_data)

    x_edges = range(extrema(x_data)..., length = n_bins + 1)
    y_edges = range(extrema(y_data)..., length = n_bins + 1)

    counts = zeros(Float64, n_bins, n_bins)
    for (x, y) in zip(x_data, y_data)
        x_bin = clamp(searchsortedfirst(x_edges, x) - 1, 1, n_bins)
        y_bin = clamp(searchsortedfirst(y_edges, y) - 1, 1, n_bins)
        counts[x_bin, y_bin] += 1
    end

    log_counts = log10.(counts .+ 1)

    # Nouveau colormap avec fond blanc, et transition de tons pastels vers des couleurs plus foncées
    custom_colormap = cgrad([:white, "#c7e9b4", "#7fcdbb", "#41b6c4", "#253494"], [0.0, 0.25, 0.5, 0.75, 1.0])

    y_max = maximum(y_data)

    with_theme(theme_latexfonts()) do
        fig = Figure(resolution = (600, 600))
        ax = Axis(
            fig[1, 1],
            title = title,
            xlabel = xlabel,
            ylabel = ylabel,
            xgridvisible = false,
            ygridvisible = false,
            xticklabelsize = 22,
            yticklabelsize = 22,
            xlabelsize = 25,
            ylabelsize = 25,
        )

        ylims!(minimum(y_data), y_max)
        hm = heatmap!(ax, x_edges, y_edges, log_counts, colormap = custom_colormap)
        Colorbar(fig[1, 2], hm, label = colorbar_label, width = 15, labelsize = 25, ticklabelsize = 22)

        if T_values !== nothing
            plot_isotherms(ax, x_data, T_values, ["CNM ($(T_values[1]) K)", "WNM ($(T_values[2]) K)"], [:blue, :red], y_max, percentage=percentage)
        end

        fig
    end
end

function plot_heatmap(data, x_values, y_values, file_path; colormap=:viridis, colorbar_label="", xlabel="x-axis (pix)", ylabel="y-axis (pix)")
    fig = Figure()
    ax = Axis(fig[1,1], xlabel=xlabel, ylabel=ylabel)
    hm = heatmap!(ax, x_values, y_values, data, colormap=colormap)
    Colorbar(fig[:, end+1], hm, label=colorbar_label)
    save(file_path, fig)
end