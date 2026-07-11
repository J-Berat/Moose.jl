from __future__ import annotations

import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

from python import moose_frontend


class MooseFrontendStateTests(unittest.TestCase):
    def parse_options(self, argv: list[str]) -> moose_frontend.FrontendOptions:
        parser = moose_frontend.build_parser()
        parsed = parser.parse_args(argv)
        return moose_frontend.validate_args(parser, parsed)

    def test_validation_normalizes_without_mutating_raw_namespace(self) -> None:
        parser = moose_frontend.build_parser()
        parsed = parser.parse_args(["--los", "Z,all", "--dry-run"])

        options = moose_frontend.validate_args(parser, parsed)

        self.assertEqual(parsed.los, ["Z,all"])
        self.assertEqual(options.los, ("z", "x", "y"))

    def test_build_julia_args_uses_validated_state(self) -> None:
        options = self.parse_options(
            [
                "--simu",
                "/data/sim-a",
                "--los",
                "x,y",
                "--faraday",
                "y",
                "--phimin",
                "-10",
                "--phimax",
                "10",
                "--dphi",
                "0.5",
                "--precision",
                "float32",
                "--tile-size",
                "64",
                "--resume",
                "safe",
                "--plan",
                "--density-kind",
                "mass_density",
                "--mean-molecular-weight",
                "1.4",
                "--hydrogen-mass-g",
                "1.6726231e-24",
                "--quiet",
            ]
        )

        self.assertEqual(
            moose_frontend.build_julia_args(options),
            [
                "--simu",
                "/data/sim-a",
                "--los",
                "x",
                "--los",
                "y",
                "--density-kind",
                "mass_density",
                "--mean-molecular-weight",
                "1.4",
                "--hydrogen-mass-g",
                "1.6726231e-24",
                "--faraday",
                "Y",
                "--phimin",
                "-10.0",
                "--phimax",
                "10.0",
                "--dphi",
                "0.5",
                "--precision",
                "float32",
                "--tile-size",
                "64",
                "--resume",
                "safe",
                "--quiet",
                "--plan",
            ],
        )

    def test_dry_run_prints_command_without_requiring_julia(self) -> None:
        options = self.parse_options(
            [
                "--julia-binary",
                "definitely-not-installed-julia",
                "--simu",
                "/data/sim-a",
                "--los",
                "z",
                "--dry-run",
            ]
        )

        stdout = io.StringIO()
        with patch.object(moose_frontend.shutil, "which", side_effect=AssertionError):
            with redirect_stdout(stdout):
                moose_frontend.run_julia_frontend(options)

        self.assertIn("Julia command:", stdout.getvalue())
        self.assertIn("definitely-not-installed-julia", stdout.getvalue())
        self.assertIn(f"--project={moose_frontend.REPO_ROOT}", stdout.getvalue())
        self.assertIn(str(moose_frontend.JULIA_ENTRYPOINT), stdout.getvalue())

    def test_dry_run_logs_composed_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "invocations.jsonl"
            options = self.parse_options(
                [
                    "--julia-binary",
                    "definitely-not-installed-julia",
                    "--los",
                    "all",
                    "--dry-run",
                    "--log-file",
                    str(log_path),
                ]
            )

            with redirect_stdout(io.StringIO()):
                moose_frontend.run_julia_frontend(options)

            log_text = log_path.read_text(encoding="utf-8")
            self.assertIn('"status": 0', log_text)
            self.assertIn('"message": "Dry run: Julia command not executed."', log_text)
            self.assertIn('"--los", "x", "--los", "y", "--los", "z"', log_text)


if __name__ == "__main__":
    unittest.main()
