"""
Prompting helpers shared by the interactive workflows.
"""

function ask_user(prompt::String, default::Float64)
    while true
        println(prompt, "(default: ", default, "): ")
        val = strip(readline())
        isempty(val) && return default

        parsed = tryparse(Float64, val)
        parsed === nothing && println("[Warning] Please enter a numeric value (e.g., 1.0) or press Enter to use the default.")
        parsed !== nothing && return parsed
    end
end

function ask_user(prompt::String, default::Int64)
    while true
        println(prompt, "(default: ", default, "): ")
        val = strip(readline())
        isempty(val) && return default

        parsed = tryparse(Int, val)
        parsed === nothing && println("[Warning] Please enter an integer value (e.g., 1 or 3) or press Enter to use the default.")
        parsed !== nothing && return parsed
    end
end

function ask_user(
    prompt::String,
    default::AbstractString;
    validate::Function = _ -> true,
    error_message::AbstractString = "[Warning] Invalid input. Please try again.",
)
    while true
        println(prompt, "(default: ", default, "): ")
        response = String(strip(readline()))
        isempty(response) && (response = String(default))
        validate(response) && return response
        !isempty(error_message) && println(error_message)
    end
end
