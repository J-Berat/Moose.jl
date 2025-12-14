"""
Centralized user-facing validation helpers for filesystem inputs.
"""

export ValidationResult, ensure_directory_access, ensure_readable_file, ensure_writable_file,
       user_error_message, validation_failure, validation_success

struct ValidationResult{T}
    value::Union{T, Nothing}
    error::Union{Nothing, String}
end

validation_success(value::T) where {T} = ValidationResult{T}(value, nothing)
validation_failure(::Type{T}, msg::AbstractString) where {T} = ValidationResult{T}(nothing, String(msg))

isvalid(result::ValidationResult) = result.error === nothing

extname(path::AbstractString) = splitext(path)[2]

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
    elseif kind == :unwritable_directory
        return "[Error] Cannot write to directory '$(path)'. Check permissions."
    elseif kind == :unwritable_file
        return "[Error] Cannot write to file '$(path)'. Check permissions."
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

function ensure_writable_file(path; expected_exts::Vector{String}=String[], must_exist::Bool=false)
    expanded = expanduser(path)

    if must_exist
        validation_error = ensure_readable_file(expanded; expected_exts=expected_exts)
        validation_error === nothing || return validation_error
    elseif !isempty(expected_exts)
        ext = lowercase(extname(expanded))
        normalized = lowercase.(expected_exts)
        !(ext in normalized) && return user_error_message(:wrong_extension, expanded; expected_exts=expected_exts)
    end

    parent_dir = dirname(expanded)
    dir_validation = ensure_directory_access(parent_dir)
    dir_validation === nothing || return dir_validation

    if isfile(expanded) && !iswritable(expanded)
        return user_error_message(:unwritable_file, expanded)
    elseif !iswritable(parent_dir)
        return user_error_message(:unwritable_directory, parent_dir)
    end

    return nothing
end

