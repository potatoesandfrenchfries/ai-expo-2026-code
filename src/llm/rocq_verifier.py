"""Rocq verifier: wraps coqc subprocess for proof compilation.

Returns a structured VerificationResult with an error taxonomy so the
proof generator can give targeted correction prompts on retry.
"""

import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Optional

# Maps coqc stderr fragments to symbolic error type names
_ERROR_MAP = {
    "Error: In environment":    "TYPE_ERROR",
    "The term":                 "TYPE_ERROR",
    "Unable to unify":          "UNIFICATION_ERROR",
    "Unsatisfied goals":        "UNRESOLVED_GOALS",
    "not applicable":           "UNRESOLVED_GOALS",
    "No such goal":             "UNRESOLVED_GOALS",
    "Universe inconsistency":   "UNIVERSE_INCONSISTENCY",
    "lra failed":               "NUMERIC_BOUNDS",
    "lia failed":               "NUMERIC_BOUNDS",
    "Cannot find":              "IMPORT_ERROR",
    "Require":                  "IMPORT_ERROR",
}


@dataclass
class VerificationResult:
    success: bool
    proof: str
    error_lines: list
    error_type: Optional[str]


class RocqVerifier:
    """Compiles a Rocq proof string with coqc and returns a VerificationResult."""

    def __init__(self, rocq_lib_path: str = "src/rocq"):
        self.rocq_lib_path = rocq_lib_path

    def verify(self, proof_text: str) -> VerificationResult:
        # Write proof to a temp file inside rocq_lib_path so -R resolves imports
        os.makedirs(self.rocq_lib_path, exist_ok=True)
        fd, tmp_path = tempfile.mkstemp(suffix=".v", dir=self.rocq_lib_path)
        try:
            with os.fdopen(fd, "w") as f:
                f.write(proof_text)

            result = subprocess.run(
                ["coqc", "-R", self.rocq_lib_path, "Chemistry", tmp_path],
                capture_output=True,
                text=True,
                timeout=30,
            )

            if result.returncode == 0:
                return VerificationResult(
                    success=True,
                    proof=proof_text,
                    error_lines=[],
                    error_type=None,
                )

            error_lines = result.stderr.strip().splitlines()
            error_type = self._classify_error(result.stderr)
            return VerificationResult(
                success=False,
                proof=proof_text,
                error_lines=error_lines,
                error_type=error_type,
            )

        except subprocess.TimeoutExpired:
            return VerificationResult(
                success=False,
                proof=proof_text,
                error_lines=["coqc timed out after 30 s"],
                error_type="TIMEOUT",
            )
        except FileNotFoundError:
            return VerificationResult(
                success=False,
                proof=proof_text,
                error_lines=["coqc not found in PATH"],
                error_type="COQC_NOT_FOUND",
            )
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

    def _classify_error(self, stderr: str) -> str:
        for fragment, label in _ERROR_MAP.items():
            if fragment in stderr:
                return label
        return "UNKNOWN_ERROR"
