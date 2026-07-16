#!/usr/bin/env julia

using Pluto

const NOTEBOOK = joinpath(@__DIR__, "MOOSE_tutorial.jl")
const OUTPUT = joinpath(@__DIR__, "MOOSE_tutorial.html")

session = Pluto.ServerSession()
notebook = Pluto.SessionActions.open(session, NOTEBOOK; run_async = false)
errors = filter(cell -> cell.errored, notebook.cells)

if !isempty(errors)
    for cell in errors
        message = get(cell.output.body, :plain_error, string(cell.output.body))
        println(stderr, cell.cell_id, " | ", first(split(message, '\n')))
    end
    error("Notebook export aborted: $(length(errors)) cell(s) failed.")
end

write(OUTPUT, Pluto.generate_html(notebook))
println("Exported $(length(notebook.cells)) cells to $OUTPUT ($(filesize(OUTPUT)) bytes).")
