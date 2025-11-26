# MOOSE

**Mock Observation Of Synchrotron Emission** is an interactive Julia toolkit for processing mock synchrotron emission from MHD simulations. It guides you through selecting simulations, configuring physical units, and running processing pipelines that compute Stokes parameters, rotation measures, and Faraday dispersion functions.

---

## Table of contents
- [Main features](#main-features)
- [Project layout](#project-layout)
- [Installation](#installation)
- [Usage](#usage)
- [Input data requirements](#input-data-requirements)
- [Outputs](#outputs)
- [Contributors](#contributors)

---

## Main features
- Compute synchrotron Stokes parameters **I**, **Q**, and **U** with optional Faraday rotation and filtering.
- Run **Rotation Measure Synthesis (RM Synthesis)** to explore Faraday depth structure.
- Generate Faraday dispersion functions, polarization angles, and derived statistics.
- Configure instrumental parameters (frequency coverage, box size, kernel filtering) interactively.
- Reuse previous answers via a saved `moose_config.json` so repeated runs only prompt for missing values.

---

## Project layout
Key source files live under `src/`:
- `MOOSE.jl`: module entrypoint that includes all domain-specific components and re-exports the main run helpers.
- `SyntheticObservations/MOOSE.jl`: interactive entrypoint that orchestrates the full workflow.
- `FileIO/ReadSimulation.jl`: helpers for loading simulation cubes and metadata.
- `Synchrotron/ProcessSynchrotron.jl`: computations for Stokes parameters and emissivity interpolation.
- `Faraday/RMSynthesis.jl`: utilities for RM synthesis and Faraday depth handling.
- `Filtering/Filter.jl` and `Utils/*.jl`: common filtering, plotting, and logging helpers.

---

## Installation
1. Install **Julia 1.10+** from [julialang.org](https://julialang.org/downloads/).
2. From the repository root, instantiate the environment to install dependencies. If you keep a personal `~/.julia/config/startup.jl`, disable it so it does not pull in packages (e.g. `JSON`) that are not part of this project:
   ```bash
   julia --startup-file=no --project -e 'using Pkg; Pkg.instantiate()'
   ```

> **Note:** Some containerized or CI environments (including the one used for automated linting here) do not ship with Julia by default. In those cases, install Julia first or run commands on a machine where Julia is available before attempting to instantiate or test the project.

---

## Usage
1. Start Julia with the project activated (again disabling any personal startup file if it imports extra packages):
   ```bash
   julia --startup-file=no --project
   ```
2. Load the interactive tool and launch it:
   ```julia
   using MOOSE
   MOOSE()
   ```
3. Follow the prompts to choose simulations, set unit conversions, select lines of sight, and enable options such as Faraday rotation or filtering. Defaults are provided for each prompt, and previous answers are reused when `moose_config.json` is present.

### Tips
- The emissivity interpolation file defaults to `Synchrotron/emissivity.dat` under your home directory; provide a different path when prompted if needed.
- Use the `help=true` keyword (`MOOSE(help=true)`) to print a detailed description of available options without running the pipeline.

### How to confirm things work
- **Smoke test the installation:** run `julia --project -e 'using MOOSE; MOOSE(help=true)'`. This precompiles the package and prints the built-in help without needing any data files.
- **Interactive end-to-end run:** follow the standard `MOOSE()` workflow described above with a real simulation directory. Success is indicated by FITS outputs next to your simulation files and a `MOOSE_summary.log` entry summarizing the run.
- **Config-driven batch run:** prepare a JSON config (for example by saving answers from a previous interactive session) and run `julia --project src/MOOSE_cli.jl /path/to/config.json --quiet`. This reuses stored parameters and will append to `MOOSE_summary.log` on completion.

---

## Recommendations
- **Start with the smoke test** before pointing to simulation data so you know the environment is healthy and dependencies are precompiled.
- **Keep `moose_config.json` under version control** (or copy it alongside each dataset) to reuse validated parameter choices and document the provenance of your outputs.
- **Name FITS files exactly as expected** (`Bx.fits`, `By.fits`, `Bz.fits`, `density.fits`, `temperature.fits`, `densityHp.fits`) to avoid interactive prompts failing on missing inputs.
- **Record the `MOOSE_summary.log` and generated FITS products together** so downstream analysis can reference both the data and the processing history.
- **Use the CLI for repeatable runs** (`julia --project src/MOOSE_cli.jl <config>.json`) and reserve the interactive session for initial exploration or parameter tuning.
- **Run on datasets stored locally** when possible; large FITS cubes streamed over networked filesystems can slow down interpolation and Faraday synthesis steps.

---

## Input data requirements
MOOSE expects simulation outputs in a directory containing the following FITS cubes with these exact filenames:
- `Bx.fits`, `By.fits`, `Bz.fits`: magnetic field components (µG recommended).
- `density.fits`: neutral hydrogen number density (cm⁻³ recommended).
- `temperature.fits`: gas temperature (K).
- `densityHp.fits`: optional electron density cube when providing `n_e` directly; otherwise it is derived from prescriptions you choose during prompts.

Rename your files before running if they use different names, for example:
```bash
mv mag_field_x.fits Bx.fits
mv n.fits density.fits
```

---

## Outputs
- Processed maps (Stokes parameters, RM maps, Faraday dispersion functions) are written alongside the simulations you choose.
- A summary log `MOOSE_summary.log` is saved in the base directory, capturing simulations processed, lines of sight, and timing information.
- Configuration choices are cached in `moose_config.json` to streamline subsequent runs.

---

## Contributors
- **Jack Berat** — Main developer
