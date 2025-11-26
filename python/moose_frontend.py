"""
Python front-end to invoke the Julia-based MOOSE CLI.

The script mirrors the flags available in `src/MOOSE_cli.jl` and forwards them
to a Julia process. It allows users who prefer Python tooling to prepare a
configuration, tweak overrides, and launch the MOOSE pipeline without having
to memorize the Julia command invocation.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List

REPO_ROOT = Path(__file__).resolve().parent.parent
JULIA_ENTRYPOINT = REPO_ROOT / "src" / "MOOSE_cli.jl"


class JuliaInvocationError(RuntimeError):
    """Raised when the Julia subprocess returns a non-zero exit code."""


def build_julia_args(parsed: argparse.Namespace) -> List[str]:
    args: List[str] = []

    config_path = parsed.config or parsed.config_path
    if config_path:
        args.append(str(Path(config_path)))

    if parsed.base_dir:
        args.extend(["--base-dir", str(parsed.base_dir)])

    for simu in parsed.simu or []:
        args.extend(["--simu", str(simu)])

    for los in parsed.los or []:
        args.extend(["--los", los.lower()])

    if parsed.interpolation:
        args.extend(["--interpolation", str(parsed.interpolation)])

    if parsed.conversionB is not None:
        args.extend(["--conversionB", str(parsed.conversionB)])

    if parsed.conversionn is not None:
        args.extend(["--conversionn", str(parsed.conversionn)])

    if parsed.conversionT is not None:
        args.extend(["--conversionT", str(parsed.conversionT)])

    if parsed.faraday:
        args.extend(["--faraday", parsed.faraday.upper()])

    if parsed.phimin is not None:
        args.extend(["--phimin", str(parsed.phimin)])

    if parsed.phimax is not None:
        args.extend(["--phimax", str(parsed.phimax)])

    if parsed.dphi is not None:
        args.extend(["--dphi", str(parsed.dphi)])

    if parsed.filtering:
        args.extend(["--filtering", parsed.filtering.upper()])

    if parsed.kernel_size is not None:
        args.extend(["--kernel-size", str(parsed.kernel_size)])

    if parsed.noise:
        args.extend(["--noise", parsed.noise.upper()])

    if parsed.snr is not None:
        args.extend(["--snr", str(parsed.snr)])

    if parsed.ne_option:
        args.extend(["--ne-option", str(parsed.ne_option)])

    if parsed.quiet:
        args.append("--quiet")

    return args


def run_julia_frontend(parsed: argparse.Namespace) -> None:
    cmd = [parsed.julia_binary, "--project", str(JULIA_ENTRYPOINT)]
    cmd.extend(build_julia_args(parsed))

    process = subprocess.run(cmd, cwd=REPO_ROOT)
    if process.returncode != 0:
        raise JuliaInvocationError(
            f"Julia front-end exited with status {process.returncode}."
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Python wrapper for the Julia MOOSE CLI. If no config is provided, "
            "the Julia tool will prompt for missing values and save them."
        )
    )

    parser.add_argument(
        "config",
        nargs="?",
        help="Optional path to a JSON config; passed as the first positional argument.",
    )
    parser.add_argument(
        "--config",
        dest="config_path",
        help="Explicit path to a JSON config file (equivalent to the positional argument).",
    )
    parser.add_argument(
        "--julia-binary",
        default="julia",
        help="Julia executable to invoke (default: julia in PATH).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress status output from the Julia pipeline.",
    )

    parser.add_argument("--base-dir", help="Base directory containing the simulations.")
    parser.add_argument(
        "--simu",
        action="append",
        help="Simulation path (may be provided multiple times).",
    )
    parser.add_argument(
        "--los",
        action="append",
        choices=["x", "y", "z", "X", "Y", "Z"],
        help="Line of sight options; may be repeated.",
    )
    parser.add_argument(
        "--interpolation",
        help="Path to an emissivity interpolation file.",
    )
    parser.add_argument("--conversionB", type=float, help="Magnetic field conversion factor.")
    parser.add_argument("--conversionn", type=float, help="Density conversion factor.")
    parser.add_argument("--conversionT", type=float, help="Temperature conversion factor.")
    parser.add_argument(
        "--faraday",
        choices=["Y", "N", "y", "n"],
        help="Enable or disable Faraday rotation.",
    )
    parser.add_argument("--phimin", type=float, help="Minimum Faraday depth.")
    parser.add_argument("--phimax", type=float, help="Maximum Faraday depth.")
    parser.add_argument("--dphi", type=float, help="Faraday depth resolution.")
    parser.add_argument(
        "--filtering",
        choices=["Y", "N", "y", "n"],
        help="Enable or disable synchrotron filtering.",
    )
    parser.add_argument("--kernel-size", dest="kernel_size", type=float, help="Filtering kernel size.")
    parser.add_argument(
        "--noise",
        choices=["Y", "N", "y", "n"],
        help="Enable or disable noise injection.",
    )
    parser.add_argument("--snr", type=float, help="Signal-to-noise ratio for noise injection.")
    parser.add_argument(
        "--ne-option",
        choices=["1", "2", "3"],
        help="Electron density prescription option (1, 2, or 3).",
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    parsed = parser.parse_args(argv)
    run_julia_frontend(parsed)


if __name__ == "__main__":
    try:
        main()
    except JuliaInvocationError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError as exc:
        missing = exc.filename if exc.filename else "julia"
        print(f"Unable to find required executable: {missing}", file=sys.stderr)
        sys.exit(1)
