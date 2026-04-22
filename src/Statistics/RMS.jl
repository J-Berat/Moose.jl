"""
    RMS(first::AbstractArray, others::AbstractArray...)

Compute the root-mean-square of one or more arrays after centering each
input around its mean. All inputs are flattened so that the statistic is
independent of array shape, and they must all contain the same number of
elements.
"""
function RMS(first::AbstractArray, others::AbstractArray...)
    vectors = [vec(A) for A in (first, others...)]

    lengths = length.(vectors)
    first_len = lengths[1]
    if any(l -> l != first_len, Iterators.drop(lengths, 1))
        throw(ArgumentError("All inputs to RMS must have the same number of elements; got lengths $(join(lengths, ", "))."))
    end

    centered_squares = ((v .- mean(v)).^2 for v in vectors)
    total = reduce(+, centered_squares)

    return sqrt(mean(total))
end
