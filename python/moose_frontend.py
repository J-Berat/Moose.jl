"""
Python front-end to invoke the Julia-based MOOSE CLI.

The script mirrors the flags available in `src/MOOSE_cli.jl` and forwards them
to a Julia process. It allows users who prefer Python tooling to prepare a
configuration, tweak overrides, and launch the MOOSE pipeline without having
to memorize the Julia command invocation.
"""

from __future__ import annotations

import argparse
import json
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable, List

REPO_ROOT = Path(__file__).resolve().parent.parent
JULIA_ENTRYPOINT = REPO_ROOT / "src" / "MOOSE_cli.jl"


class JuliaInvocationError(RuntimeError):
    """Raised when the Julia subprocess returns a non-zero exit code."""


def resolve_config_path(parsed: argparse.Namespace) -> Path | None:
    """Return the selected config path, favoring the positional argument."""

    config_value = parsed.config or parsed.config_path
    return Path(config_value) if config_value else None


def build_julia_args(parsed: argparse.Namespace, config_path: Path | None) -> List[str]:
    args: List[str] = []

    if config_path:
        args.append(str(Path(config_path)))

    if parsed.base_dir:
        args.extend(["--base-dir", str(parsed.base_dir)])

    for simu in parsed.simu or []:
        args.extend(["--simu", str(simu)])

    for los in parsed.los or []:
        args.extend(["--los", los])

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

    if parsed.zeta is not None:
        args.extend(["--zeta", str(parsed.zeta)])

    if parsed.Geff is not None:
        args.extend(["--Geff", str(parsed.Geff)])

    if parsed.phiPAH is not None:
        args.extend(["--phiPAH", str(parsed.phiPAH)])

    if parsed.XC is not None:
        args.extend(["--XC", str(parsed.XC)])

    if parsed.quiet:
        args.append("--quiet")

    return args


def format_command(command: Iterable[str]) -> str:
    """Return a shell-safe string representation of the command."""

    return shlex.join(list(command))


def _normalize_los_values(values: List[str] | None, parser: argparse.ArgumentParser) -> list[str]:
    """Return a normalized list of LOS choices, mirroring the Julia CLI semantics."""

    if not values:
        return []

    normalized: list[str] = []
    seen: set[str] = set()
    for raw_value in values:
        parts = [value.strip().lower() for value in raw_value.split(",") if value.strip()]
        for value in parts:
            if value == "all":
                for expanded in ("x", "y", "z"):
                    if expanded not in seen:
                        normalized.append(expanded)
                        seen.add(expanded)
                continue

            if value not in {"x", "y", "z"}:
                parser.error("--los must be x, y, z, or 'all'.")

            if value not in seen:
                normalized.append(value)
                seen.add(value)

    return normalized


def validate_args(parser: argparse.ArgumentParser, parsed: argparse.Namespace) -> Path | None:
    """Validate CLI arguments and return the resolved config path."""

    if parsed.config and parsed.config_path:
        parser.error("Provide either a positional config path or --config, not both.")

    parsed.los = _normalize_los_values(parsed.los, parser)
    config_path = resolve_config_path(parsed)
    if config_path and not config_path.exists():
        parser.error(f"Config file not found: {config_path}")

    if parsed.filtering and parsed.filtering.upper() == "Y" and parsed.kernel_size is None:
        parser.error("--kernel-size is required when enabling filtering.")

    if parsed.noise and parsed.noise.upper() == "Y" and parsed.snr is None:
        parser.error("--snr is required when enabling noise injection.")

    if parsed.faraday and parsed.faraday.upper() == "Y":
        missing_faraday = [
            flag
            for flag, value in {
                "--phimin": parsed.phimin,
                "--phimax": parsed.phimax,
                "--dphi": parsed.dphi,
            }.items()
            if value is None
        ]
        if missing_faraday:
            parser.error(
                "Faraday rotation enabled; please provide values for "
                + ", ".join(missing_faraday)
            )

    return config_path


def _log_invocation(log_file: Path, command: Iterable[str], status: int, message: str | None) -> None:
    """Append a JSON log entry describing the invocation outcome."""

    log_file.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "command": list(command),
        "command_string": format_command(command),
        "status": status,
        "message": message,
    }
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry))
        handle.write("\n")


def run_julia_frontend(parsed: argparse.Namespace) -> None:
    config_path = resolve_config_path(parsed)
    cmd = [parsed.julia_binary, "--project", str(JULIA_ENTRYPOINT)]
    cmd.extend(build_julia_args(parsed, config_path))

    julia_path = shutil.which(parsed.julia_binary)
    if julia_path is None:
        raise FileNotFoundError(parsed.julia_binary)

    command_string = format_command(cmd)
    if parsed.print_command:
        print(f"Julia command: {command_string}")

    if parsed.dry_run:
        message = "Dry run: Julia command not executed."
        if parsed.log_file:
            _log_invocation(Path(parsed.log_file), cmd, status=0, message=message)
        return

    process = subprocess.run(cmd, cwd=REPO_ROOT)

    failure_message = None
    if process.returncode != 0:
        failure_message = (
            "Julia front-end exited with status "
            f"{process.returncode}. Command: {command_string}"
        )

    if parsed.log_file:
        _log_invocation(Path(parsed.log_file), cmd, process.returncode, failure_message)

    if failure_message:
        raise JuliaInvocationError(failure_message)


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
    parser.add_argument(
        "--print-command",
        action="store_true",
        help="Display the composed Julia command before running it.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the Julia command and exit without invoking Julia.",
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
        help=(
            "Line of sight options; may be repeated. Accepts x, y, z, 'all', "
            "or comma-separated combinations such as x,y,z."
        ),
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
    parser.add_argument("--zeta", type=float, help="Wolfire model zeta constant.")
    parser.add_argument("--Geff", type=float, help="Wolfire model Geff constant.")
    parser.add_argument("--phiPAH", type=float, help="Wolfire model phiPAH constant.")
    parser.add_argument("--XC", type=float, help="Wolfire model XC constant.")
    parser.add_argument(
        "--log-file",
        help="Optional path to write a JSONL log entry for each invocation.",
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    parsed = parser.parse_args(argv)
    validate_args(parser, parsed)
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
