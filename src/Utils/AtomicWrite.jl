"""
    atomic_write_path(path, writer)

Write a file through a temporary sibling and atomically move it into place.
The temporary file is created in the destination directory so the final `mv`
stays on the same filesystem.
"""
function atomic_write_path(path::AbstractString, writer::Function)
    final_path = abspath(expanduser(String(path)))
    parent_dir = dirname(final_path)
    mkpath(parent_dir)

    _, ext = splitext(final_path)
    tmp_path = tempname(parent_dir) * ext

    try
        writer(tmp_path)
        mv(tmp_path, final_path; force = true)
    catch
        rm(tmp_path; force = true)
        rethrow()
    end

    return final_path
end

atomic_write_path(writer::Function, path::AbstractString) = atomic_write_path(path, writer)

function atomic_write_text(path::AbstractString, content::AbstractString)
    return atomic_write_path(path) do tmp_path
        open(tmp_path, "w") do io
            write(io, content)
        end
    end
end
