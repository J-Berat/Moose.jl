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
2. From the repository root, instantiate the environment to install dependencies:
   ```bash
   julia --project -e 'using Pkg; Pkg.instantiate()'
   ```

---

## Usage
1. Start Julia with the project activated:
   ```bash
   julia --project
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
