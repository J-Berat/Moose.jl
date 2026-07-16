#!/usr/bin/env julia

using Pluto

const NOTEBOOK = joinpath(@__DIR__, "MOOSE_tutorial.jl")
const TWO_SCREENS = "Synchrotron background + two Faraday screens"

function failed_cells(notebook)
    filter(cell -> cell.errored, notebook.cells)
end

function print_failures(notebook, scenario)
    for cell in failed_cells(notebook)
        message = get(cell.output.body, :plain_error, string(cell.output.body))
        println(stderr, "[$scenario] ", cell.cell_id, " | ", first(split(message, '\n')))
    end
end

session = Pluto.ServerSession()
notebook = Pluto.SessionActions.open(session, NOTEBOOK; run_async = false)
print_failures(notebook, "one screen")
isempty(failed_cells(notebook)) || error("Default one-screen notebook failed.")

notebook.bonds[:mock_faraday_case] = Dict("value" => TWO_SCREENS)
Pluto.set_bond_values_reactive(
    ; session, notebook, bound_sym_names = [:mock_faraday_case],
    is_first_values = [false], run_async = false,
)
print_failures(notebook, "two screens")
isempty(failed_cells(notebook)) || error("Two-screen notebook failed.")

validation_cell = only(filter(
    cell -> string(cell.cell_id) == "00000000-0000-0000-0000-000000001180",
    notebook.cells,
))
occursin("❌", validation_cell.output.body) && error("Two-screen validation reported a failed check.")

println("Notebook checks passed for one-screen and two-screen mock scenarios.")
