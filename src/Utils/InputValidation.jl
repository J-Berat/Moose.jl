"""
Centralized user-facing validation helpers for filesystem inputs.
"""

export ensure_directory_access, ensure_readable_file, user_error_message

function user_error_message(kind::Symbol, path; expected_exts::Vector{String}=String[], reason::Union{Nothing, AbstractString}=nothing)
    if kind == :missing_directory
        return "[Error] The directory '$(path)' does not exist."
    elseif kind == :unreadable_directory
        return "[Error] Cannot read from directory '$(path)'. Check permissions."
    elseif kind == :missing_file
        return "[Error] The file '$(path)' does not exist."
    elseif kind == :not_regular_file
        return "[Error] The path '$(path)' is not a regular file."
    elseif kind == :unreadable_file
        return "[Error] Cannot read the file '$(path)'. Check permissions."
    elseif kind == :wrong_extension
        expected = join(expected_exts, " or ")
        return "[Error] The file '$(path)' does not match the expected extension ($(expected))."
    elseif kind == :corrupted_file
        detail = reason === nothing ? "" : " Details: $(reason)."
        return "[Error] The file '$(path)' could not be read and may be corrupted." * detail
    else
        return "[Error] Invalid input at '$(path)'."
    end
end

function ensure_directory_access(path)
    if !isdir(path)
        return user_error_message(:missing_directory, path)
    elseif !isreadable(path)
        return user_error_message(:unreadable_directory, path)
    end

    return nothing
end

function ensure_readable_file(path; expected_exts::Vector{String}=String[])
    if !ispath(path)
        return user_error_message(:missing_file, path)
    elseif !isfile(path)
        return user_error_message(:not_regular_file, path)
    end

    if !isempty(expected_exts)
        ext = lowercase(extname(path))
        normalized = lowercase.(expected_exts)
        !(ext in normalized) && return user_error_message(:wrong_extension, path; expected_exts=expected_exts)
    end

    if !isreadable(path)
        return user_error_message(:unreadable_file, path)
    end

    return nothing
end

