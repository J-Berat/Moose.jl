"""
Common error helpers for user-facing entrypoints.
"""

export MooseError, cli_error, config_error, throw_cli_error, throw_config_error

struct MooseError <: Exception
    code::Symbol
    message::String
    exit_code::Int

    function MooseError(code::Symbol, message::AbstractString; exit_code::Int=1)
        return new(code, String(message), exit_code)
    end
end

Base.showerror(io::IO, err::MooseError) = print(io, err.message)

cli_error(message::AbstractString; code::Symbol=:cli_invalid_argument, exit_code::Int=2) =
    MooseError(code, "[CLI Error] " * String(message); exit_code=exit_code)

config_error(message::AbstractString; code::Symbol=:config_invalid, exit_code::Int=3) =
    MooseError(code, "[Config Error] " * String(message); exit_code=exit_code)

throw_cli_error(message::AbstractString; code::Symbol=:cli_invalid_argument, exit_code::Int=2) =
    throw(cli_error(message; code=code, exit_code=exit_code))

throw_config_error(message::AbstractString; code::Symbol=:config_invalid, exit_code::Int=3) =
    throw(config_error(message; code=code, exit_code=exit_code))
