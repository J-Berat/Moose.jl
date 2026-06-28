# MOOSE

**Mock Observation Of Synchrotron Emission** is an interactive Julia toolkit for processing mock synchrotron emission from MHD simulations. It guides you through selecting simulations, configuring physical units, and running processing pipelines that compute Stokes parameters, rotation measures, and Faraday dispersion functions. The goal is to make it straightforward to turn raw simulation cubes into reproducible FITS products that mirror common radio-observational analyses.

---

## Table of contents
- [Main features](#main-features)
- [Project layout](#project-layout)
- [Quickstart](#quickstart)
- [Installation](#installation)
- [Docker](#docker)
- [Usage](#usage)
- [Public API](#public-api)
- [Recommendations](#recommendations)
- [Configuration file schemas](#configuration-file-schemas)
- [Input data requirements](#input-data-requirements)
- [Outputs](#outputs)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Citation](#citation)
- [Contributors](#contributors)

---

## Main features
- Compute synchrotron Stokes parameters **I**, **Q**, and **U** with optional Faraday rotation and filtering.
- Run **Rotation Measure Synthesis (RM Synthesis)** to explore Faraday depth structure.
- Generate Faraday dispersion functions, polarization angles, and derived statistics.
- Configure instrumental parameters (frequency coverage, box size, interferometric Fourier filtering) interactively.
- Run either interactively (`run_moose`) or non-interactively from JSON config (`src/MOOSE_cli.jl` / `MOOSE_from_config`).

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

## Quickstart
If you already have Julia 1.10+ installed, the fastest way to confirm the tool works is:

```bash
julia --startup-file=no --project -e 'using Pkg; Pkg.instantiate(); using MOOSE; run_moose(help=true)'
```

This installs dependencies, precompiles the package, and prints the interactive help without requiring any input data. Once that succeeds, move on to preparing your simulation directory and running the standard workflow in [Usage](#usage).

---

## Installation
1. Install **Julia 1.10+** from [julialang.org](https://julialang.org/downloads/).
2. From the repository root, instantiate the environment to install dependencies. If you keep a personal `~/.julia/config/startup.jl`, disable it so it does not pull in packages (e.g. `JSON`) that are not part of this project:
   ```bash
   julia --startup-file=no --project -e 'using Pkg; Pkg.instantiate()'
   ```

   This project tracks `Manifest.toml` for reproducible dependency resolution. Running `Pkg.instantiate()` will install the pinned dependency set for your platform.

> **Note:** Some containerized or CI environments (including the one used for automated linting here) do not ship with Julia by default. In those cases, install Julia first or run commands on a machine where Julia is available before attempting to instantiate or test the project.

---

## Docker
Build the Julia container from the repository root:

```bash
docker build -t moose-julia .
```

The image uses Julia 1.12.6 to match the checked-in `Manifest.toml`.

The default command prints the built-in help:

```bash
docker run --rm moose-julia
```

For a config-driven run, mount the directory that contains your simulation data and config file, then call the Julia CLI inside the container:

```bash
docker run --rm \
  -v /path/to/data:/data \
  moose-julia \
  julia --startup-file=no --project=/app /app/src/MOOSE_cli.jl /data/config.json --quiet
```

The image also includes `python3`, so the Python front-end can be used if preferred:

```bash
docker run --rm \
  -v /path/to/data:/data \
  moose-julia \
  python3 /app/python/moose_frontend.py --config /data/config.json --quiet
```

---

## Usage
1. Start Julia with the project activated (again disabling any personal startup file if it imports extra packages):
   ```bash
   julia --startup-file=no --project
   ```
2. Load the interactive tool and launch it:
   ```julia
   using MOOSE
   run_moose()
   ```
3. Follow the prompts to choose simulations, set unit conversions, select lines of sight, and enable options such as Faraday rotation or filtering.

### Interactive modes
- `run_moose()` defaults to `reset_config=true`, which starts from fresh prompts.
- To reload a previous JSON config and reuse saved answers, run:
  ```julia
  run_moose(reset_config=false)
  ```

### Tips
- The emissivity interpolation file default is `~/emissivity.dat`; provide a different path when prompted if needed.
- Use the `help=true` keyword (`run_moose(help=true)`) to print a detailed description of available options without running the pipeline.
- Keep the Julia session running while iterating on parameters; repeated calls to `run_moose` reuse precompiled code and are noticeably faster than restarting Julia each time.

### Julia CLI (non-interactive)
Use the Julia CLI entrypoint when you already have a config file or want scripted runs:

```bash
julia --startup-file=no --project src/MOOSE_cli.jl /path/to/config.json --quiet
```

You can also provide or override values directly with flags such as `--base-dir`, `--simu`, `--los`, `--interpolation`, `--faraday`, `--filtering`, `--noise`, and `--ne-option`.

Important behavior:
- Without `--write-back`, CLI overrides are merged in-memory. The summary log records the config that was read, the effective in-memory config label, and that no config was saved.
- With `--write-back`, the provided config JSON is overwritten with merged values before execution.
- This CLI is non-interactive: missing required values (for example `base_dir`, simulations, or interpolation path) produce an explicit error instead of prompts.
- `--filtering Y --kernel-size <L>` applies a hard Fourier-domain 0/1 interferometric mask. The value `<L>` is the largest retained spatial scale in pixels, matching the filtering convention used in the Depolarization instrumental pipeline.

### Python front-end
Prefer Python tooling? Use the lightweight wrapper in `python/moose_frontend.py`, which forwards familiar CLI flags to the Julia entrypoint:

```bash
python python/moose_frontend.py --simu /data/simulation --los z --quiet
```

The wrapper accepts the same options documented for `src/MOOSE_cli.jl` (for example, `--conversionB`, `--filtering`, `--ne-option`, `--rng-seed`, and positional or `--config` paths). The `--julia-binary` flag lets you point to a non-default Julia executable when needed.

By default, supplying a config file does **not** overwrite that file; overrides are applied in memory. Add `--write-back` if you explicitly want to persist overrides into the provided config JSON.

For quick validation, `--print-command` shows the fully composed Julia invocation before running it, and `--dry-run` prints the command and exits without launching Julia. You can also add `--log-file /path/to/invocations.jsonl` to append one JSONL entry per run.

Validation rules in the Python wrapper:
- `--faraday Y` requires `--phimin`, `--phimax`, and `--dphi`.
- `--filtering Y` requires `--kernel-size`.
- `--noise Y` requires `--snr`.
- `--write-back` requires a config path (`--config` or positional).

### How to test locally
Run the test suite to validate the installation and catch regressions before processing large datasets:

```bash
julia --startup-file=no --project -e 'using Pkg; Pkg.test()'
```

### How to confirm things work
- **Smoke test the installation:** run `julia --startup-file=no --project -e 'using MOOSE; run_moose(help=true)'`. This precompiles the package and prints the built-in help without needing any data files.
- **Interactive end-to-end run:** follow the standard `run_moose()` workflow described above with a real simulation directory. Success is indicated by FITS outputs next to your simulation files and a `MOOSE_summary.log` entry summarizing the run.
- **Config-driven batch run:** prepare a JSON config (for example by saving answers from a previous interactive session) and run `julia --project src/MOOSE_cli.jl /path/to/config.json --quiet`. This reuses stored parameters and will append to `MOOSE_summary.log` on completion. Pass `--write-back` only when you want CLI overrides persisted into that config file.
- **Use the provided template:** copy `config/default_config.json`, update paths/constants (`base_dir`, `simulations`, conversion factors, frequencies, Faraday/noise/filter flags), then run it with the CLI. The config loader now validates required fields and fails fast with explicit errors if a value is invalid.

---

## Public API
The stable Julia API is the set of names exported by `using MOOSE`:

- `run_moose` for the interactive workflow.
- `MOOSE_from_config` for JSON-driven batch runs.
- `MooseError`, `cli_error`, and `config_error` for user-facing failures.
- Faraday tomography: `RMClean`, `RMCleanHealpix`, `RMCleanAuto`, `RMCleanResult`, `rmsf_diagnostics`, `RMSFDiagnostics`, and `write_rmsf`.
- HEALPix helpers: `HealpixStack`, `HealpixRMResult`, `RMSynthesisHealpix`, `RMSynthesisAuto`, `healpix_map`, `healpix_maps_from_stack`, `detect_fits_grid`, `is_healpix_fits`, `is_image_fits`, `read_fits_grid`, `read_fits_grid_stack`, `read_healpix_map`, `read_healpix_stack`, `write_healpix_map`, `write_healpix_stack`, and `write_healpix_rm_result`.

Other qualified names such as `MOOSE.RMS`, `MOOSE.Pnu`, or `MOOSE.buildHeader3D` are internal implementation details. They are tested for regression coverage, but they are not yet promised as stable external API.

---

## Recommendations
- **Start with the smoke test** before pointing to simulation data so you know the environment is healthy and dependencies are precompiled.
- **Keep `moose_config.json` under version control** (or copy it alongside each dataset) to reuse validated parameter choices and document the provenance of your outputs.
- **Name FITS files exactly as expected** (`Bx.fits`, `By.fits`, `Bz.fits`, `density.fits`, `temperature.fits`, `densityHp.fits`) to avoid interactive prompts failing on missing inputs.
- **Record the `MOOSE_summary.log` and generated FITS products together** so downstream analysis can reference both the data and the processing history.
- **Use the CLI for repeatable runs** (`julia --project src/MOOSE_cli.jl <config>.json`) and reserve the interactive session for initial exploration or parameter tuning.
- **Run on datasets stored locally** when possible; large FITS cubes streamed over networked filesystems can slow down interpolation and Faraday synthesis steps.

### User pre-run checklist
1. Run a quick pre-check before a full run:
   ```bash
   julia --startup-file=no --project -e 'using MOOSE; run_moose(help=true)'
   ```
2. Keep config files read-only by default: pass a config JSON without `--write-back`, and only use `--write-back` when you explicitly want to persist CLI overrides.
3. Validate command composition before execution (Python wrapper):
   ```bash
   python3 python/moose_frontend.py --config cfg.json --print-command --dry-run
   ```
4. Keep one folder per dataset, storing together: `config.json`, `MOOSE_summary.log`, and the generated FITS outputs.
5. Verify required input filenames before running: `Bx.fits`, `By.fits`, `Bz.fits`, `density.fits`, `temperature.fits`, `densityHp.fits`.

---

## Configuration file schemas
`MOOSE_from_config` supports two JSON styles:

1. Flat keys (legacy/interactive style):
   - `base_dir`, `simulations` or `chosen_simu`, `chosen_LOS`
   - `conversionB`, `conversionn`, `conversionT`
   - `FaradayRotation`, `phimin`, `phimax`, `dphi`
   - `do_rm_clean`, `rm_clean_gain`, `rm_clean_niter`, `rm_clean_threshold`
   - `responseSynchrotron`, `kernel_size_synchrotron` (largest Fourier scale retained by the interferometric 0/1 mask, in pixels)
   - `add_noise`, `SNR_nu`
   - `interpolation_file_path`
   - `ne_option`, `IonizationFraction`
   - `BoxLength_pc`, `BoxLength_pix`
   - `nustart`, `nuend`, `dnu`

2. Nested keys (template/frontend style), as in `config/default_config.json`:
   - `freq.start|end|step`
   - `box.size_pc|npix` for cubic boxes, or `box.x|y|z` plus `box.npix` for LOS-dependent physical sizes
   - `faraday.enabled|phimin|phimax|dphi`
   - `rm_clean.enabled|gain|niter|threshold`
   - `emissivity.path`
   - `ne.mode|ion_fraction`
   - `rng_seed` for reproducible noise injection

Notes:
- `base_dir` is required.
- Simulations are required (`simulations` or `chosen_simu`).
- Emissivity path is required (`interpolation_file_path` or `emissivity.path`).
- Relative simulation paths and emissivity paths are resolved against `base_dir`.
- Config frequency values (`nustart`, `nuend`, `dnu`, or `freq.start|end|step`) are in MHz.
- FITS spectral axes are written in Hz (`CUNIT3 = "Hz"`), after a single MHz-to-Hz conversion at write time.
- Faraday depth values are in `rad/m^2`.
- RM-CLEAN requires Faraday rotation to be enabled and a uniformly spaced Faraday-depth grid.
- Box lengths are in parsec and pixel counts are dimensionless.

---

## Input data requirements
MOOSE expects simulation outputs in a directory containing the following FITS cubes with these exact filenames:
- `Bx.fits`, `By.fits`, `Bz.fits`: magnetic field components (µG).
- `density.fits`: neutral hydrogen number density (cm⁻³).
- `temperature.fits`: gas temperature (K).
- `densityHp.fits`: optional electron density cube when providing `n_e` directly; otherwise it is derived from prescriptions you choose during prompts.

Rename your files before running if they use different names, for example:
```bash
mv mag_field_x.fits Bx.fits
mv n.fits density.fits
```

### HEALPix maps

MOOSE also exposes helpers for HEALPix maps through [Healpix.jl](https://github.com/JuliaAstro/Healpix.jl). This is useful when Q/U observations or mock products are already stored as one HEALPix FITS map per frequency:

```julia
using MOOSE

q = read_healpix_stack(q_files)  # Npix x Nfreq
u = read_healpix_stack(u_files)

nu_hz = [120e6, 121e6, 122e6]
phi = -50.0:0.5:50.0

result = RMSynthesisHealpix(q, u, nu_hz, phi)
write_healpix_rm_result("healpix_rm", result; prefix="lofar")
```

All maps in a stack must have the same `NSIDE`, ordering (`RING` or `NESTED`), and frequency order. You can also pass plain `Npix x Nfreq` matrices with `RMSynthesisHealpix(Q, U, nu_hz, phi; order=:ring)`.

MOOSE can also detect the FITS grid type automatically. A regular FITS image/cube is routed to the normal cube routines, while a HEALPix FITS binary table (`PIXTYPE=HEALPIX`) is routed through the HEALPix path:

```julia
detect_fits_grid("Qnu.fits")          # :image or :healpix

# Regular cube paths, or vectors of HEALPix map paths ordered like nu_hz.
result = RMSynthesisAuto(q_input, u_input, nu_hz, phi)
clean = RMCleanAuto(q_input, u_input, nu_hz, phi; gain=0.1, niter=1000)
```

For HEALPix frequency stacks, pass one FITS file per frequency, for example `q_input = ["Q_0001.fits", "Q_0002.fits", ...]`. For cube data, pass the single Q and U cube paths.

The full `MOOSE_from_config` / `run_moose_processing` pipeline also detects HEALPix physical inputs automatically. Keep the same field names as the cube workflow (`Bx`, `By`, `Bz`, `density`, `temperature`, and optional `densityHp`). Each `*.fits` may be a single HEALPix map or a HEALPix cube/stack stored in one binary table, either as multiple map columns or as one vector column with one row per HEALPix pixel:

```text
simulation/
  Bx.fits
  By.fits
  Bz.fits
  density.fits
  temperature.fits
```

You can also use one directory per field containing ordered LOS shells:

```text
simulation/
  Bx/shell_0001.fits
  Bx/shell_0002.fits
  By/shell_0001.fits
  ...
```

All physical fields in one simulation must use the same grid type. HEALPix inputs preserve `NSIDE` and ordering in the outputs, and MOOSE writes integrated maps plus frequency/Faraday stacks as HEALPix FITS products. The current `responseSynchrotron=Y` Fourier filter is cartesian only; HEALPix runs should keep it disabled until a spherical-harmonic filter is added.

### RMSF diagnostics and RM-CLEAN

RM synthesis returns a *dirty* Faraday dispersion function (FDF): the true Faraday spectrum convolved with the Rotation Measure Spread Function (RMSF). Two helpers characterise and deconvolve it.

`rmsf_diagnostics(nu_hz, phi)` returns an `RMSFDiagnostics` object with the complex RMSF on a symmetric lag grid plus the standard resolution metrics — the measured Faraday resolution `fwhm` (δφ), the analytic `fwhm_theoretical = 2√3/Δλ²`, the maximum recoverable Faraday depth `phi_max`, and the largest sensitive Faraday scale `max_scale`. These are written automatically as `RMSF.fits` (with the metrics in the header) whenever the pipeline runs with Faraday rotation enabled.

`RMClean` runs RM synthesis and then deconvolves the FDF with an RM-CLEAN loop (Heald 2009), restoring the clean components with a Gaussian beam matched to the RMSF main lobe:

```julia
using MOOSE

# Q, U with frequency as the last axis (1D/2D/3D); frequencies in Hz.
result = RMClean(Q, U, nu_hz, phi; gain=0.1, niter=2000, threshold=1e-3)

result.cleanFDF    # |restored FDF|
result.model       # clean-component model
result.residual    # residual FDF
result.rmsf.fwhm   # restoring-beam FWHM (rad/m^2)
```

For HEALPix inputs, `RMCleanHealpix(Q, U, nu_hz, phi; order=:ring)` returns a `HealpixRMResult` of the restored FDF that can be written directly with `write_healpix_rm_result`. `phi` must be uniformly spaced; its step sets the RMSF sampling.

In config-driven runs, enable RM-CLEAN with:

```json
"rm_clean": {
  "enabled": true,
  "gain": 0.1,
  "niter": 1000,
  "threshold": 0.0
}
```

When enabled, the pipeline writes `cleanFDF.fits`, `realCleanFDF.fits`, `imagCleanFDF.fits`, and `residualFDF.fits` next to the dirty FDF products. HEALPix runs write one map per Faraday-depth slice for each of those quantities.

---

## Outputs
- Processed maps (Stokes parameters, RM maps, Faraday dispersion functions) are written alongside the simulations you choose.
- Each synchrotron output directory also includes polarization diagnostic PNGs for the brightest `Pnumax` sightline: `polarization_angle_vs_lambda2.png`, `fractional_polarization_vs_lambda2.png`, and `stokes_qu_diagram.png`.
- When Faraday rotation is enabled, an `RMSF.fits` file is written next to the FDF products, holding the complex RMSF (`|R|`, `Re R`, `Im R`) and the resolution metrics (`RMSFFWHM`, `RMSFTHEO`, `PHIMAX`, `MAXSCALE`) in its header.
- When RM-CLEAN is enabled, restored and residual FDF products are written as `cleanFDF`, `realCleanFDF`, `imagCleanFDF`, and `residualFDF`.
- HEALPix RM-synthesis outputs are written as one standard HEALPix FITS map per Faraday-depth slice.
- A summary log `MOOSE_summary.log` is saved in the base directory, capturing simulations processed, lines of sight, and timing information.
- Configuration choices are cached in `moose_config.json` to streamline subsequent runs.
- FITS headers include reproducibility metadata such as MOOSE version, git hash, config hash, config provenance, line of sight, run options, and key unit conventions.

---

## Troubleshooting
- **Julia not found:** verify `julia --version` works in your shell or pass `--julia-binary` to the Python wrapper.
- **Pkg errors about missing system dependencies:** rerun `Pkg.instantiate()` with `--startup-file=no` to avoid loading packages from a global startup file, and ensure internet access is available for downloads.
- **Cannot locate FITS cubes:** double-check filenames and working directory; the interactive prompts will echo the expected paths when they are missing.
- **CLI exits with missing required fields:** ensure your JSON defines `base_dir`, at least one simulation (`simulations`/`chosen_simu`), and an emissivity path (`interpolation_file_path`/`emissivity.path`).
- **Slow repeated runs:** keep `moose_config.json` alongside each dataset to skip prompts and reuse validated parameters.

---

## License
MOOSE is distributed under the MIT License. See [LICENSE](LICENSE).

---

## Citation
If you use MOOSE, please cite the associated paper: [2026A&A...708A.245B](https://ui.adsabs.harvard.edu/abs/2026A%26A...708A.245B/abstract).

---

## Contributors
- **Jack Berat** — Main developer
