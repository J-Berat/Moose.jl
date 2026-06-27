"""
Prompting helpers shared by the interactive workflows.
"""

# --- internal formatting helpers -------------------------------------------

const _PROMPT_MARK = "▸"

# Strip a trailing colon / whitespace so we can append " [default]: " cleanly,
# regardless of how the caller phrased the prompt.
_clean_prompt(prompt::AbstractString) = rstrip(rstrip(prompt), ':') |> rstrip

function _print_prompt(prompt::AbstractString, default)
    printstyled("  ", _PROMPT_MARK, " "; color = :cyan, bold = true)
    print(_clean_prompt(prompt))
    printstyled(" [", default, "]"; color = :light_cyan)
    print(": ")
    flush(stdout)
end

function _warn(msg::AbstractString)
    printstyled("    ↳ ", msg, "\n"; color = :yellow)
    flush(stdout)
end

# --- public API ------------------------------------------------------------

function ask_user(prompt::String, default::Float64)
    while true
        _print_prompt(prompt, default)
        val = strip(readline())
        isempty(val) && return default

        parsed = tryparse(Float64, val)
        parsed === nothing && _warn("Please enter a numeric value (e.g., 1.0) or press Enter to use the default.")
        parsed !== nothing && return parsed
    end
end

function ask_user(prompt::String, default::Int64)
    while true
        _print_prompt(prompt, default)
        val = strip(readline())
        isempty(val) && return default

        parsed = tryparse(Int, val)
        parsed === nothing && _warn("Please enter an integer value (e.g., 1 or 3) or press Enter to use the default.")
        parsed !== nothing && return parsed
    end
end

function ask_user(
    prompt::String,
    default::AbstractString;
    validate::Function = _ -> true,
    error_message::AbstractString = "Invalid input. Please try again.",
)
    while true
        _print_prompt(prompt, default)
        response = String(strip(readline()))
        isempty(response) && (response = String(default))
        validate(response) && return response
        !isempty(error_message) && _warn(error_message)
    end
end
