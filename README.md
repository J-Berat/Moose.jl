# Moose

**Moose** (*Mock Observation Of Synchrotron Emission*) is an interactive Julia toolkit for processing mock synchrotron emission from MHD simulations. It guides you through selecting simulations, configuring physical units, and running processing pipelines that compute Stokes parameters, rotation measures, and Faraday dispersion functions. The goal is to make it straightforward to turn raw simulation cubes into reproducible FITS products that mirror common radio-observational analyses.

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
- Fit physical Faraday models directly to `q(λ²)`/`u(λ²)` spectra (**QU fitting**: external screen, Burn slab, external/internal Faraday dispersion) with AIC/BIC model selection, per spectrum or per pixel.
- Generate Faraday dispersion functions, polarization angles, and derived statistics.
- Turbulence diagnostics: polarization gradient maps **|∇P|** (Gaensler et al. 2011) and isotropic **structure functions** of RM or polarization-angle maps (π-ambiguity aware).
- Produce per-pixel **spectral index maps** (`alpha.fits`, `alpha_err.fits`) from the multi-frequency intensity cube via a log-log least-squares fit.
- Configure instrumental parameters (frequency coverage, box size, interferometric Fourier filtering) interactively.
- Run either interactively (`run_moose`) or non-interactively from JSON config (`src/MOOSE_cli.jl` / `MOOSE_from_config`).
- Generate a self-contained **demo dataset with analytically known results** (`make_demo_data`) to validate an installation end to end.

---

## Project layout
Key source files live under `src/`:
- `Moose.jl`: module entrypoint that includes all domain-specific components and re-exports the main run helpers.
- `SyntheticObservations/Moose.jl`: interactive entrypoint that orchestrates the full workflow.
- `FileIO/ReadSimulation.jl`: helpers for loading simulation cubes and metadata.
- `Synchrotron/ProcessSynchrotron.jl`: computations for Stokes parameters and emissivity interpolation.
- `Faraday/RMSynthesis.jl`: utilities for RM synthesis and Faraday depth handling.
- `Filtering/Filter.jl` and `Utils/*.jl`: common filtering, plotting, and logging helpers.

---

## Quickstart
If you already have Julia 1.10+ installed, the fastest way to confirm the tool works is:

```bash
julia --startup-file=no --project -e 'using Pkg; Pkg.instantiate(); using Moose; run_moose(help=true)'
```

This installs dependencies, precompiles the package, and prints the interactive help without requiring any input data. Once that succeeds, move on to preparing your simulation directory and running the standard workflow in [Usage](#usage).

### Demo dataset with known results
To run the full pipeline without any real simulation data, generate the built-in demo dataset and process it:

```julia
using Moose
demo = make_demo_data("moose_demo")
MOOSE_from_config(demo.config_path; quiet = true)
```

The demo is a Faraday screen in front of a uniform power-law emitter, built so that every output is **analytically known**: `demo.expected` holds the exact values (`rm` for `RMmap.fits`, `alpha` for `alpha.fits`, `Tnu` per channel, `qnu_over_tnu`/`unu_over_tnu`, the polarization fraction, `intne`, `intBLOS`, and the FDF peak position). The same values are written to `moose_demo/expected_results.json` so you can compare the generated FITS products against them by hand. If they match, your installation is processing data correctly end to end.

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
   using Moose
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
- **Smoke test the installation:** run `julia --startup-file=no --project -e 'using Moose; run_moose(help=true)'`. This precompiles the package and prints the built-in help without needing any data files.
- **Interactive end-to-end run:** follow the standard `run_moose()` workflow described above with a real simulation directory. Success is indicated by FITS outputs next to your simulation files and a `MOOSE_summary.log` entry summarizing the run.
- **Config-driven batch run:** prepare a JSON config (for example by saving answers from a previous interactive session) and run `julia --project src/MOOSE_cli.jl /path/to/config.json --quiet`. This reuses stored parameters and will append to `MOOSE_summary.log` on completion. Pass `--write-back` only when you want CLI overrides persisted into that config file.
- **Use the provided template:** copy `config/default_config.json`, update paths/constants (`base_dir`, `simulations`, conversion factors, frequencies, Faraday/noise/filter flags), then run it with the CLI. The config loader now validates required fields and fails fast with explicit errors if a value is invalid.

---

## Public API
The stable Julia API is the set of names exported by `using Moose`:

- `run_moose` for the interactive workflow.
- `MOOSE_from_config` for JSON-driven batch runs.
- `MooseError`, `cli_error`, and `config_error` for user-facing failures.
- Faraday tomography: `RMClean`, `RMCleanHealpix`, `RMCleanAuto`, `RMCleanResult`, `rmsf_diagnostics`, `RMSFDiagnostics`, and `write_rmsf`.
- QU fitting: `QUFit` (single spectrum), `QUFitCompare` (fit all models, rank by BIC), `QUFitCube` (per-pixel parameter/uncertainty/χ²red maps, NaN-masked pixels skipped), `QUFitResult`, `qu_model`, and `QU_FIT_MODELS` (`:screen`, `:burn_slab`, `:external_dispersion`, `:internal_dispersion`). Inputs follow the `RMSynthesis` conventions (frequency in Hz on the last axis); optional `sigma_q`/`sigma_u` uncertainties weight the fit and calibrate `chi2_red` and parameter errors.
- Spectral index: `spectral_index_map` for per-pixel log-log power-law fits of intensity cubes.
- Turbulence diagnostics: `polarization_gradient_map` (per-map or per-channel `|∇P|`, optional `normalized=|∇P|/|P|`, `pixel_size` scaling, NaN-masked pixels preserved) and `structure_function` / `StructureFunctionResult` (Monte-Carlo pair sampling with logarithmic separation bins; pass `angle=true` for polarization-angle maps so differences are wrapped modulo π; seedable `rng` for reproducibility).
- Demo dataset: `make_demo_data` to generate a synthetic quickstart dataset with analytically known results.
- HEALPix helpers: `HealpixStack`, `HealpixRMResult`, `RMSynthesisHealpix`, `RMSynthesisAuto`, `healpix_map`, `healpix_maps_from_stack`, `detect_fits_grid`, `is_healpix_fits`, `is_image_fits`, `read_fits_grid`, `read_fits_grid_stack`, `read_healpix_map`, `read_healpix_stack`, `read_healpix_cube`, `write_healpix_map`, `write_healpix_stack`, `write_healpix_cube`, `write_healpix_rm_result`, `healpix_udgrade`, `healpix_reorder`, `healpix_smooth`, and `HEALPIX_UNSEEN`. Masked pixels (UNSEEN sentinel) are converted to `NaN` when reading stacks and back to the sentinel on write (opt out with `unseen_to_nan`/`nan_to_unseen`); masked pixels are skipped by RM synthesis and RM-CLEAN. `COORDSYS` metadata is preserved end-to-end, and `write_healpix_stack`/`write_healpix_rm_result` accept `format=:cube` to write a single multi-slice FITS file (with a `COORDS` extension) instead of one file per slice. Simulation fields with mismatched NSIDE or ordering are automatically conformed (with a warning) to the `Bx` reference grid via `healpix_udgrade`/`healpix_reorder`; `healpix_smooth` applies a Gaussian beam in spherical-harmonic space (`fwhm_arcmin`/`fwhm_deg`/`fwhm_rad`, NaN-aware masked smoothing). `tile_size` also works with HEALPix inputs (single-column HEALPix tables): the sky is processed in bands of HEALPix pixel rows and 3D products are streamed as single-file HEALPix cubes.
- **HEALPix vector-field convention (important):** for HEALPix simulations the line of sight is the local radial direction of each pixel, so only `LOS = "z"` is accepted (`"x"`/`"y"` raise a config error). The `Bx`/`By`/`Bz` (and `Vx`/`Vy`/`Vz`) FITS inputs must contain the **per-pixel tangent-basis components** — `Bx = B·e_θ` (colatitude), `By = B·e_φ` (azimuth), `Bz = B·e_r` (radial/LOS) — not global cartesian components. The intrinsic polarization angle `ψ = atan(B2, B1) + π/2` is then measured in the local `(e_θ, e_φ)` basis. If your simulation stores global cartesian vectors, project them onto each pixel's tangent basis before writing the FITS files.

Other qualified names such as `Moose.RMS`, `Moose.Pnu`, or `Moose.buildHeader3D` are internal implementation details. They are tested for regression coverage, but they are not yet promised as stable external API.

---

## Recommendations
- **Start with the smoke test** before pointing to simulation data so you know the environment is healthy and dependencies are precompiled.
- **Keep `moose_config.json` under version control** (or copy it alongside each dataset) to reuse validated parameter choices and document the provenance of your outputs.
- **Use the default field names or configure `field_sources`** for inputs named differently from `Bx`, `By`, `Bz`, `density`, `temperature`, and optional `densityHp`.
- **Record the `MOOSE_summary.log` and generated FITS products together** so downstream analysis can reference both the data and the processing history.
- **Use the CLI for repeatable runs** (`julia --project src/MOOSE_cli.jl <config>.json`) and reserve the interactive session for initial exploration or parameter tuning.
- **Run on datasets stored locally** when possible; large FITS cubes streamed over networked filesystems can slow down interpolation and Faraday synthesis steps.

### User pre-run checklist
1. Run a quick pre-check before a full run:
   ```bash
   julia --startup-file=no --project -e 'using Moose; run_moose(help=true)'
   ```
2. Keep config files read-only by default: pass a config JSON without `--write-back`, and only use `--write-back` when you explicitly want to persist CLI overrides.
3. Validate command composition before execution (Python wrapper):
   ```bash
   python3 python/moose_frontend.py --config cfg.json --print-command --dry-run
   ```
4. Keep one folder per dataset, storing together: `config.json`, `MOOSE_summary.log`, and the generated FITS outputs.
5. Verify required input field sources before running: either the default names (`Bx.fits`, `By.fits`, `Bz.fits`, `density.fits`, `temperature.fits`, optional `densityHp.fits`) or the matching `field_sources` entries in your config.

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
   - `field_sources.Bx|By|Bz|density|temperature|densityHp` to map canonical field names to custom file names or HDF5 datasets
   - `physical_mask.T_min|T_max|n_min|n_max` to exclude cells outside selected temperature/density ranges
   - `density_kind`, `mean_molecular_weight`, and `hydrogen_mass_g` to interpret `density` as either number density or mass density

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
- Field source paths are resolved relative to each simulation directory unless absolute. For HDF5, use `{"file": "simulation.h5", "dataset": "group/name"}` when the dataset is not named after the canonical field.
- Physical mask thresholds are applied after unit conversions: `T_*` is in K and `n_*` is the MOOSE number density `n` in `cm^-3`.
- With `"density_kind": "number_density"` (default), `conversionn` converts the input `density` values directly to `cm^-3`. With `"density_kind": "mass_density"`, `conversionn` converts the input `density` values to `g cm^-3`, then MOOSE computes `nH = rho / (mean_molecular_weight * hydrogen_mass_g)`.

### Performance options (opt-in)
Two optional keys control the memory footprint of large runs; both default to the historical behaviour when omitted:

- `"precision": "float32"` processes and stores the cubes in single precision, halving the steady-state RAM and the size of the FITS products. Per-pixel accumulations still run in double precision, so the loss of accuracy is negligible for mock-observation work (relative errors ~1e-7). The default is `"float64"`. The FITS headers record the choice in the `PRECIS` keyword. Not supported with HEALPix inputs. RM-CLEAN products are not reduced.
- `"tile_size": N` processes the sky plane in bands of `N` rows: input cubes are read from FITS one band at a time and the 3D products (`Qnu`, `Unu`, `Tnu`, `Pnu`, `polfrac`, `ne`, FDF cubes) are streamed to disk band by band, so cubes much larger than the available RAM can be processed. The per-pixel math is identical to a plain run, so results match exactly. Recorded in the `TILESIZE` header keyword. Limitations: incompatible with interferometric filtering (`responseSynchrotron`, needs the full sky plane in Fourier space), noise injection (`add_noise`, the per-channel σ derives from the full-map rms), and RM-CLEAN; the polarization diagnostic plots are skipped in tiled runs. Both options combine for cartesian grids (`"precision": "float32"` + `"tile_size": N`).

The equivalent CLI flags are `--precision float32` and `--tile-size N` (Julia CLI and Python front-end).

Use `--plan` with either CLI to validate all required simulation fields and their dimensions without loading full cubes or producing outputs. The preflight report lists frequency/Faraday channel counts and per-LOS estimates for peak working RAM, FITS data volume, and cell-channel workload. Estimates intentionally exclude allocator overhead, plots, FITS headers, FFT workspace, and RM-CLEAN iteration cost.

- `"resume": "safe"` writes an atomic `.moose-complete.json` manifest after each simulation/LOS finishes. A later run skips that unit only when its processing configuration, input file size/timestamps, and declared FITS outputs still match. The default is `"off"`; use `--resume safe` from either CLI. Safe resume is rejected when noise injection is enabled because skipping work would otherwise alter the shared random-number sequence.
- `"outputs": ["integrated", "rm"]` restricts processing and writes to selected product groups. Available groups are `integrated` (electron-density and LOS summary maps), `stokes` (Q/U/T/P cubes and polarization fractions), `rm` (RM map), `fdf` (Faraday dispersion products and RMSF), `spectral_index` (alpha maps), and `diagnostics` (polarization plots). The default `["all"]` preserves the complete legacy run. Derived groups compute required intermediates without writing unrequested prerequisite products. Use `--outputs rm,fdf` from either CLI. `rm`/`fdf` require Faraday rotation, RM-CLEAN requires `fdf`, and selective groups are currently unavailable with tiled processing.

---

## Input data requirements
MOOSE expects simulation outputs in a directory containing the following regular cubes as FITS or HDF5 files:
- `Bx.fits`, `By.fits`, `Bz.fits` (or `.h5`/`.hdf5`): magnetic field components (µG).
- `density.fits` (or `.h5`/`.hdf5`): neutral hydrogen number density (cm⁻³).
- `temperature.fits` (or `.h5`/`.hdf5`): gas temperature (K).
- `densityHp.fits` (or `.h5`/`.hdf5`): optional electron density cube when providing `n_e` directly; otherwise it is derived from prescriptions you choose during prompts.

**Warning:** `density` means number density `n` in the MOOSE equations by default, not mass density `rho`. If your simulation stores mass density, set `density_kind` explicitly:

```json
"density_kind": "mass_density",
"mean_molecular_weight": 1.4,
"hydrogen_mass_g": 1.6726231e-24
```

In this mode, `conversionn` must convert the on-disk mass density to `g cm^-3`; MOOSE then computes `nH = rho / (mean_molecular_weight * hydrogen_mass_g)`. `densityHp`, when provided for `ne_option = 3`, remains an electron number-density cube.

To exclude cells outside a physical phase or trusted range, add a `physical_mask` block. Masked cells contribute zero to synchrotron emission, Faraday rotation, and LOS-integrated quantities:

```json
"physical_mask": {
  "T_min": 100.0,
  "T_max": 1000000.0,
  "n_min": 0.001,
  "n_max": 10.0
}
```

Aliases such as `temperature_min`, `temperature_max`, `density_min`, `density_max`, `nH_min`, and `nH_max` are also accepted in config files.

For one-file-per-field HDF5 inputs, each file may contain a single numeric dataset, or a dataset named after the field (for example `Bx` inside `Bx.h5`). You can also store all fields in one shared `.h5`/`.hdf5` file; MOOSE will look for datasets named `Bx`, `By`, `Bz`, `density`, `temperature`, and the optional density fields, including inside groups such as `/fields/Bx`.

If your files use different names, add a `field_sources` mapping to the config instead of renaming them:
```json
"field_sources": {
  "Bx": "mag_field_x.fits",
  "By": "mag_field_y.fits",
  "Bz": "mag_field_z.fits",
  "density": "nH.fits",
  "temperature": "thermal.fits",
  "densityHp": "electron_density.fits"
}
```

For shared HDF5 files with arbitrary dataset names:
```json
"field_sources": {
  "Bx": {"file": "simulation.h5", "dataset": "fields/mag_x"},
  "By": {"file": "simulation.h5", "dataset": "fields/mag_y"},
  "Bz": {"file": "simulation.h5", "dataset": "fields/mag_z"},
  "density": {"file": "simulation.h5", "dataset": "gas/nH"},
  "temperature": {"file": "simulation.h5", "dataset": "gas/temp"}
}
```

### HEALPix maps

MOOSE also exposes helpers for HEALPix maps through [Healpix.jl](https://github.com/JuliaAstro/Healpix.jl). This is useful when Q/U observations or mock products are already stored as one HEALPix FITS map per frequency:

```julia
using Moose

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
using Moose

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
- Each synchrotron output directory includes a spectral index map `alpha.fits` and its 1σ uncertainty `alpha_err.fits`, from a per-pixel log-log least-squares fit of the brightness-temperature cube `Tnu`. The map uses the flux-density convention `S_ν ∝ ν^α` (`α = β_T + 2`, recorded in the `ALPHADEF` header keyword); the uncertainty is `NaN` when only two frequency channels are available (no residual degrees of freedom).
- Each synchrotron output directory includes the polarization fraction cube `polfrac.fits` (`Pnu/Tnu`) and the map `polfracmax.fits`, the maximum finite polarization fraction over the frequency axis. Pixels with non-positive or non-finite `Tnu` are written as `NaN`.
- Each synchrotron output directory also includes polarization diagnostics for the brightest `Pnumax` sightline: individual PNGs (`polarization_angle_vs_lambda2.png`, `fractional_polarization_vs_lambda2.png`, `stokes_qu_diagram.png`) plus a publication-ready composite as `polarization_diagnostics.png` and `polarization_diagnostics.pdf`.
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
Moose is distributed under the MIT License. See [LICENSE](LICENSE).

---

## Citation
If you use Moose/MOOSE, please cite the associated paper: [2026A&A...708A.245B](https://ui.adsabs.harvard.edu/abs/2026A%26A...708A.245B/abstract).

---

## Contributors
- **Jack Berat** — Main developer
