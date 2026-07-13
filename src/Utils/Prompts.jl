"""
Prompting helpers shared by the interactive workflows.
"""

# --- internal formatting helpers -------------------------------------------

const _PROMPT_MARK = "▸"

"""Print a SHINE-style heading for a group of interactive prompts."""
function prompt_section(title::AbstractString)
    println()
    printstyled("╭─ ", title, "\n"; color = :light_cyan, bold = true)
    flush(stdout)
end

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

# Styled, homogeneous feedback used across the interactive workflows.
function warn_user(msg::AbstractString)
    printstyled("    ↳ ", msg, "\n"; color = :yellow)
    flush(stdout)
end

function error_user(msg::AbstractString)
    printstyled("    ↳ ", msg, "\n"; color = :light_red, bold = true)
    flush(stdout)
end

# Predicate for yes/no prompts: accepts only Y or N (case-insensitive),
# matching the downstream `uppercase(answer) == "Y"` checks.
is_yes_no(answer) = uppercase(strip(String(answer))) in ("Y", "N")

# --- public API ------------------------------------------------------------

function ask_user(prompt::String, default::Float64)
    while true
        _print_prompt(prompt, default)
        val = strip(readline())
        isempty(val) && return default

        parsed = tryparse(Float64, val)
        parsed === nothing && warn_user("Please enter a numeric value (e.g., 1.0) or press Enter to use the default.")
        parsed !== nothing && return parsed
    end
end

function ask_user(prompt::String, default::Int64)
    while true
        _print_prompt(prompt, default)
        val = strip(readline())
        isempty(val) && return default

        parsed = tryparse(Int, val)
        parsed === nothing && warn_user("Please enter an integer value (e.g., 1 or 3) or press Enter to use the default.")
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
        !isempty(error_message) && warn_user(error_message)
    end
end
