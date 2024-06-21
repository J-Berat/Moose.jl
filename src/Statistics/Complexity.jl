function complexity(y; x=range(1, length(y)), m=moments(y,x=x))
    
    # y : array of which the complexity is calculated on
    # x : array of indices for each value in y
    # m : moments 0, 1 and 2 of y
    #
    # MAMD 2023-11-07

    dv = abs(x[2]-x[1])

    Amp = m[1]*dv/(sqrt(2π)*m[3])
    G = gauss(x, Amp, m[2], m[3])
    delta = y .- G
    result = std(delta) / std(G)

    return result
    
end

function MapComplexity(cube, x)
    nx, ny = size(cube,1), size(cube,2)
    complexity_map = zeros(nx, ny)
    
    for i in 1:nx
        for j in 1:ny
            spectre = cube[i, j, :]
            complexity_map[i, j] = complexity(spectre, x=x, m=moments(spectre,x=x))
        end
    end

    return complexity_map
end