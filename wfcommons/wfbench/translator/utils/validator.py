"""
Skill-driven validation for LLM-translated workflow code.

Reads validation rules from the target skill's `## Validation` section and
runs them against generated code. Classifies failures as truncation (retry
with same prompt) vs. semantic errors (retry with feedback), and raises
`TranslationValidationError` when retries are exhausted.
"""

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from ..skills.loader import Skill


class TranslationValidationError(Exception):
    """Raised when generated code fails validation after all retries."""

    def __init__(self, message: str, code: str, errors: list[str]):
        super().__init__(message)
        self.code = code
        self.errors = errors


# Substrings that indicate the LLM's output was truncated mid-generation
# (typically hitting max_tokens). These are best handled with a silent retry.
_TRUNCATION_MARKERS = (
    "unterminated string literal",
    "unterminated triple-quoted",
    "EOF while scanning",
    "EOF in multi-line",
    "unexpected EOF",
)


class ValidationResult:
    """Outcome of running validation against generated code."""

    def __init__(self):
        self.errors: list[str] = []
        self.is_truncation: bool = False

    @property
    def passed(self) -> bool:
        return not self.errors

    def add_error(self, message: str, truncation: bool = False):
        self.errors.append(message)
        if truncation:
            self.is_truncation = True


def validate(code: str, skill: Skill) -> ValidationResult:
    """
    Run all skill-defined checks against the generated code.

    :param code: Generated code from the LLM.
    :param skill: The target-system skill with validation rules.
    :return: ValidationResult with errors (empty if passed).
    """
    result = ValidationResult()

    _check_required_elements(code, skill, result)
    _run_syntax_check(code, skill, result)

    return result


def _check_required_elements(code: str, skill: Skill, result: ValidationResult):
    """Verify each required element (substring) appears in the code."""
    for element in skill.required_elements:
        if element not in code:
            result.add_error(f"Missing required element: '{element}'")


def _run_syntax_check(code: str, skill: Skill, result: ValidationResult):
    """Run the skill's syntax_check_command on a temp file containing the code."""
    cmd_template = skill.syntax_check_command
    if not cmd_template:
        return

    tool = cmd_template.split()[0]
    if not shutil.which(tool) and tool != "python":
        # Tool isn't installed — skip rather than fail. This keeps validation
        # useful even when the target WMS CLI isn't available locally.
        return

    suffix = _suffix_for_skill(skill)
    with tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False) as f:
        f.write(code)
        tmp_path = f.name

    try:
        cmd = cmd_template.replace("{file}", tmp_path)
        proc = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=30)
        if proc.returncode != 0:
            err = (proc.stderr + proc.stdout).strip()
            truncation = any(m in err for m in _TRUNCATION_MARKERS)
            result.add_error(
                f"Syntax check failed ({cmd_template.split()[0]}): {err[:500]}",
                truncation=truncation,
            )
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _suffix_for_skill(skill: Skill) -> str:
    """Pick a file suffix for the temp file used in syntax checks."""
    name = skill.filepath.stem
    mapping = {
        "nextflow": ".nf",
        "cwl": ".cwl",
        "bash": ".sh",
        "airflow": ".py",
        "radical_pilot": ".py",
        "pegasus": ".py",
        "parsl": ".py",
        "dask": ".py",
        "pycompss": ".py",
        "taskvine": ".py",
        "swift_t": ".swift",
    }
    return mapping.get(name, ".txt")


def format_retry_feedback(errors: list[str]) -> str:
    """Compose a feedback block to append to a retry prompt."""
    lines = ["=== VALIDATION FEEDBACK FROM PREVIOUS ATTEMPT ===",
             "Your previous output failed validation with the following errors:"]
    for e in errors:
        lines.append(f"  - {e}")
    lines.append("Please produce a corrected version that fixes these issues.")
    return "\n".join(lines)
