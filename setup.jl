#!/usr/bin/env julia
"""
setup.jl

Simple helper to prepare the MOOSE environment. It activates the project in
this repository, installs dependencies, and precompiles the code. Pass
`--test` to also run the package tests after setup.

Usage examples (disable personal startup files that import extra packages):
  julia --startup-file=no setup.jl
  julia --startup-file=no setup.jl --test
"""

using Pkg

function main(; run_tests::Bool=false)
    project_dir = @__DIR__

    println("Activating project at $(abspath(project_dir))…")
    Pkg.activate(project_dir)

    println("Instantiating dependencies (this may take a while if this is the first run)…")
    Pkg.instantiate()

    println("Precompiling project…")
    Pkg.precompile()

    if run_tests
        println("Running test suite…")
        Pkg.test()
    else
        println("Skipping tests (pass --test to run them).")
    end

    println("Moose is ready. You can now launch it with `using Moose; run_moose()`.")
end

main(run_tests = "--test" in ARGS)
