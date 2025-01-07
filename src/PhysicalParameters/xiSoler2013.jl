function compute_xi(A_parallel, A_perp)
    if (A_parallel + A_perp) == 0
        error("The sum of A_parallel and A_perp cannot be zero.")
    end
    return (A_parallel - A_perp) / (A_parallel + A_perp)
end

function calculate_areas(data, bin_edges)
    # Create the histogram
    hist = fit(Histogram, data, bin_edges)
    
    # Extract bin edges and counts
    edges = hist.edges[1]
    counts = hist.weights

    # Define ranges for A_parallel and A_perp
    A_parallel_range = (cos_phi -> 0.75 < cos_phi <= 1)
    A_perp_range = (cos_phi -> 0 < cos_phi <= 0.25)
    
    # Integrate over the ranges
    A_parallel = sum(counts[i] for i in 1:length(counts) if A_parallel_range(edges[i]))
    A_perp = sum(counts[i] for i in 1:length(counts) if A_perp_range(edges[i]))
    
    return A_parallel, A_perp
end