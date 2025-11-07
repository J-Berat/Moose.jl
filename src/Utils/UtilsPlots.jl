function plot_isotherms!(ax, x_data, T_values, labels, colors, y_max; percentage = 0.35, show_isotherms = true)
    if show_isotherms && T_values !== nothing
        sorted_T_values = sort(T_values)
        for (i, T) in enumerate(sorted_T_values)
            n_range = range(minimum(x_data), stop = maximum(x_data), length = 100)
            P_line = log10.(10 .^ n_range .* T)

            lines!(ax, n_range, P_line,
                   linestyle = :dash,
                   linewidth = 3,
                   color = colors[i],
                   label = "T = $T K")

            pos_index = Int(percentage * length(n_range))
            x_pos = Float32(n_range[pos_index + 1])
            y_pos = Float32(P_line[pos_index - 3])
            delta_x = n_range[pos_index + 7] - n_range[pos_index - 7]
            delta_y = P_line[pos_index + 15] - P_line[pos_index - 15]
            angle = atan(delta_y, delta_x) .- deg2rad(28)

            text!(
                Point(x_pos, y_pos),
                text = labels[i],
                color = colors[i],
                align = (:center, :center),
                rotation = angle,
                fontsize = 12
            )
        end
    end
end

function phase_diagram(x_data, y_data;
    T_values = nothing,
    n_bins = 500,
    title = "",
    xlabel = "",
    ylabel = "",
    colorbar_label = "",
    percentage = 0.35,
    apply_log = true,
    show_isotherms = true,
    show_Tequ = true,
    colormap = cgrad([:white, "#c7e9b4", "#7fcdbb", "#41b6c4", "#253494"], [0.0, 0.25, 0.5, 0.75, 1.0]),
    savepath = nothing)

    with_theme(theme_latexfonts()) do

        fig = Figure(backgroundcolor = :white, resolution = (800, 600))
    
        # Safety check for log scaling
        if apply_log && (any(x_data .<= 0) || any(y_data .<= 0))
            error("x_data and y_data must be strictly positive before log transformation.")
        end
    
        # Log-transform if needed
        if apply_log
            x_data = log10.(x_data)
            y_data = log10.(y_data)
        end
    
        # Axis
        ax = Axis(fig[1, 1],
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            xticklabelsize = 30,
            yticklabelsize = 30,
            xlabelsize = 30,
            ylabelsize = 30,
            xgridvisible = false,
            ygridvisible = false
        )
    
        # Bin edges
        x_min, x_max = extrema(x_data)
        y_min, y_max = extrema(y_data)
    
        x_edges = range(x_min, x_max, length = n_bins + 1)
        y_edges = range(y_min, y_max, length = n_bins + 1)
    
        # Histogram faster version
        counts = zeros(Float64, n_bins, n_bins)
    
        x_idx = clamp.(searchsortedfirst.(Ref(x_edges), x_data) .- 1, 1, n_bins)
        y_idx = clamp.(searchsortedfirst.(Ref(y_edges), y_data) .- 1, 1, n_bins)
    
        for i in eachindex(x_idx)
            counts[x_idx[i], y_idx[i]] += 1
        end
    
        log_counts = log10.(counts .+ 1)
    
        # Plot the heatmap
        hm = heatmap!(ax, x_edges, y_edges, log_counts; colormap = colormap)
    
        # Fix limits
        xlims!(ax, x_min, x_max)
        ylims!(ax, y_min, y_max)
    
        # Colorbar placed on fig[1,2]
        if colorbar_label != ""
            cb = Colorbar(fig[1, 2], hm, label = colorbar_label, ticklabelsize = 30, labelsize = 30)
        else
            cb = Colorbar(fig[1, 2], hm, ticklabelsize = 30, labelsize = 30)
        end
    
        # Make colorbar closer
        colgap!(fig.layout, 5)
    
        # Optional isotherms
        if show_isotherms && T_values !== nothing
            labels = ["", ""]
            colors = [:blue, :red]
            y_max = maximum(y_data)
            plot_isotherms!(ax, x_data, T_values, labels, colors, y_max; percentage = percentage, show_isotherms = show_isotherms)
        end
    
        # Optional Tequilibrium line
        if show_Tequ
            nequ = logindgen(1000, 0.01, 1000)
            Tequ = @. tequilibrium(nequ)
            Pequilibre = nequ .* Tequ
    
            lines!(ax, log10.(nequ), log10.(Pequilibre),
                linestyle = :dot,
                linewidth = 4,
                color = :black,
                label = "Thermal equilibrium"
            )
        end
    
        # Minor ticks settings for log scales
        if apply_log
            ax.xminorgridvisible = false
            ax.yminorgridvisible = false
            ax.xminorticksvisible = true
            ax.yminorticksvisible = true
        end
    
        # Legend
        leg = axislegend(ax, position = :rb, labelsize = 30)
        leg.framevisible = false
    
        # Save if needed
        if savepath !== nothing
            save(savepath, fig)
        end
    
        return fig
    end
end

function phase_diagram!(ax::Axis, x_data, y_data;
    T_values = nothing,
    n_bins = 500,
    percentage = 0.35,
    apply_log = true,
    show_isotherms = true,
    show_Tequ = true,
    colormap = cgrad([:white, "#c7e9b4", "#7fcdbb", "#41b6c4", "#253494"], [0.0, 0.25, 0.5, 0.75, 1.0])
)
    if apply_log && (any(x_data .<= 0) || any(y_data .<= 0))
        error("x_data and y_data must be strictly positive before log transformation.")
    end

    if apply_log
        x_data = log10.(x_data)
        y_data = log10.(y_data)
    end

    x_min, x_max = extrema(x_data)
    y_min, y_max = extrema(y_data)

    x_edges = range(x_min, x_max, length = n_bins + 1)
    y_edges = range(y_min, y_max, length = n_bins + 1)

    counts = zeros(Float64, n_bins, n_bins)
    x_idx = clamp.(searchsortedfirst.(Ref(x_edges), x_data) .- 1, 1, n_bins)
    y_idx = clamp.(searchsortedfirst.(Ref(y_edges), y_data) .- 1, 1, n_bins)

    for i in eachindex(x_idx)
        counts[x_idx[i], y_idx[i]] += 1
    end

    log_counts = log10.(counts .+ 1)

    hm = heatmap!(ax, x_edges, y_edges, log_counts; colormap = colormap)

    xlims!(ax, x_min, x_max)
    ylims!(ax, y_min, y_max)

    if show_isotherms && T_values !== nothing
        plot_isotherms!(ax, x_data, T_values, ["200 K", "2000 K"], [:blue, :red], y_max;
                        percentage = percentage, show_isotherms = show_isotherms)
    end

    if show_Tequ
        nequ = logindgen(1000, 0.01, 1000)
        Tequ = @. tequilibrium(nequ)
        Pequilibre = nequ .* Tequ
        lines!(ax, log10.(nequ), log10.(Pequilibre), linestyle = :dot, linewidth = 3, color = :black)
    end

    return hm
end

function hist2D(x_data, y_data;
    n_bins = 100,
    apply_log = false,
    title = "",
    xlabel = "",
    ylabel = "",
    colorbar_label = "",
    colormap = cgrad([:white, "#c7e9b4", "#7fcdbb", "#41b6c4", "#253494"], [0.0, 0.25, 0.5, 0.75, 1.0]),
    savepath = nothing)

    fig = Figure(backgroundcolor = :white, resolution = (800, 600))

    if apply_log && (any(x_data .<= 0) || any(y_data .<= 0))
        error("x_data and y_data must be strictly positive before log transformation.")
    end

    if apply_log
        x_data = log10.(x_data)
        y_data = log10.(y_data)
    end

    x_min, x_max = extrema(x_data)
    y_min, y_max = extrema(y_data)

    x_edges = range(x_min, x_max, length = n_bins + 1)
    y_edges = range(y_min, y_max, length = n_bins + 1)

    counts = zeros(Float64, n_bins, n_bins)
    x_idx = clamp.(searchsortedfirst.(Ref(x_edges), x_data) .- 1, 1, n_bins)
    y_idx = clamp.(searchsortedfirst.(Ref(y_edges), y_data) .- 1, 1, n_bins)

    for i in eachindex(x_idx)
        counts[x_idx[i], y_idx[i]] += 1
    end

    values = apply_log ? log10.(counts .+ 1) : counts

    ax = Axis(fig[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticklabelsize = 16,
        yticklabelsize = 16,
        xlabelsize = 20,
        ylabelsize = 20
    )

    hm = heatmap!(ax, x_edges, y_edges, values; colormap = colormap)

    xlims!(ax, x_min, x_max)
    ylims!(ax, y_min, y_max)

    if colorbar_label != ""
        Colorbar(fig[1, 2], hm, label = colorbar_label, ticklabelsize = 16)
    else
        Colorbar(fig[1, 2], hm, ticklabelsize = 16)
    end

    colgap!(fig.layout, 5)

    if savepath !== nothing
        save(savepath, fig)
    end

    return fig
end