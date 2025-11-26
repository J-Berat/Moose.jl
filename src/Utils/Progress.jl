"""
Simple progress bar used when iterating over simulations.
"""

function print_progress(progress::Int, total::Int)
    @assert total > 0 "Total must be greater than 0."
    @assert 0 <= progress <= total "Progress must be between 0 and total."

    bar_width = 50
    progress_ratio = progress / total
    filled_length = Int(round(bar_width * progress_ratio))
    empty_length = bar_width - filled_length

    green = "\u001b[42m"
    reset = "\u001b[0m"
    gray = "\u001b[47m"

    filled_bar = green * repeat(" ", filled_length) * reset
    empty_bar = gray * repeat(" ", empty_length) * reset

    percentage = Int(round(100 * progress_ratio))

    print("\rProgress: |$filled_bar$empty_bar| $progress/$total ($percentage%)")
    flush(stdout)

    if progress == total
        println()
    end
end
