"""
    RMS(first::AbstractArray, others::AbstractArray...)

Compute the root-mean-square of one or more arrays after centering each
input around its mean. All inputs are flattened so that the statistic is
independent of array shape.
"""
function RMS(first::AbstractArray, others::AbstractArray...)
    centered_squares = ((vec(A) .- mean(vec(A))).^2 for A in (first, others...))
    total = reduce(+, centered_squares)
    return sqrt(mean(total))
end

