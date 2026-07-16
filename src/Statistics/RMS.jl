"""
    RMS(first::AbstractArray, others::AbstractArray...)

Compute the root-mean-square of one or more arrays after centering each
input around its mean. All inputs are flattened so that the statistic is
independent of array shape, and they must all contain the same number of
elements. Non-finite observations are ignored using a common mask across
all inputs. `NaN` is returned when no complete finite observation remains.
"""
function RMS(first::AbstractArray, others::AbstractArray...)
    vectors = [vec(A) for A in (first, others...)]

    lengths = length.(vectors)
    first_len = lengths[1]
    if any(l -> l != first_len, Iterators.drop(lengths, 1))
        throw(ArgumentError("All inputs to RMS must have the same number of elements; got lengths $(join(lengths, ", "))."))
    end

    valid = trues(first_len)
    for vector in vectors
        valid .&= isfinite.(vector)
    end
    any(valid) || return NaN

    finite_vectors = (vector[valid] for vector in vectors)
    centered_squares = ((v .- mean(v)).^2 for v in finite_vectors)
    total = reduce(+, centered_squares)

    return sqrt(mean(total))
end
