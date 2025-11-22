"""
    RMS(X)
    RMS(X, Y)
    RMS(X, Y, Z)

Compute the root-mean-square of one, two, or three arrays after centering
each input around its mean. All inputs are flattened so that the statistic
is independent of array shape.
"""
function RMS(X)
    x = vec(X)
    return sqrt(mean((x .- mean(x)).^2))
end

function RMS(X, Y)
    squares = (vec(X) .- mean(vec(X))).^2 .+ (vec(Y) .- mean(vec(Y))).^2
    return sqrt(mean(squares))
end

function RMS(X, Y, Z)
    squares = (vec(X) .- mean(vec(X))).^2 .+ (vec(Y) .- mean(vec(Y))).^2 .+
              (vec(Z) .- mean(vec(Z))).^2
    return sqrt(mean(squares))
end
