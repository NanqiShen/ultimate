#!/usr/bin/env python3
"""Assert a pytest run really passed by parsing its summary line.

A shell exit code is not acceptable proof of a green test run: a pytest usage
error (for example an unknown ``--timeout`` flag) exits non-zero but was once
swallowed by a wrapper that reported success. This script therefore reads the
captured output and requires:

* a pytest summary line matching ``N passed``,
* ``N`` greater than zero,
* no ``failed``, ``error``, ``errors`` entries,
* no ``usage:`` / ``unrecognized arguments`` diagnostic anywhere.

Exit code is 0 only when every assertion holds.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

SUMMARY_RE = re.compile(
    r"^(?=.*?\b(?P<passed>\d+) passed\b)?"
    r"(?=.*?\b(?P<failed>\d+) failed\b)?"
    r"(?=.*?\b(?P<error>\d+) errors?\b)?"
    r"(?=.*?\b(?P<skipped>\d+) skipped\b)?"
    r".*?\bin (?P<duration>[\d.]+)s\b"
)


def parse_summary(text: str) -> dict[str, int] | None:
    """Return counts from the pytest summary line, or None when absent."""
    for line in reversed(text.splitlines()):
        stripped = line.strip()
        if not stripped or " in " not in stripped:
            continue
        if not re.search(r"\b(passed|failed|error|errors|skipped|no tests ran)\b", stripped):
            continue
        if "no tests ran" in stripped:
            return {"passed": 0, "failed": 0, "error": 0, "skipped": 0, "found": 1}
        match = SUMMARY_RE.match(stripped)
        if match is None:
            continue
        groups = match.groupdict()
        return {
            "passed": int(groups["passed"] or 0),
            "failed": int(groups["failed"] or 0),
            "error": int(groups["error"] or 0),
            "skipped": int(groups["skipped"] or 0),
            "found": 1,
        }
    return None


def check(text: str) -> list[str]:
    problems: list[str] = []
    lowered = text.lower()
    if "unrecognized arguments" in lowered:
        problems.append("pytest reported unrecognized arguments (usage error), not a test result")
    if re.search(r"^usage: python", text, flags=re.MULTILINE):
        problems.append("pytest printed a usage error; no tests were executed")
    summary = parse_summary(text)
    if summary is None:
        problems.append("no pytest summary line found in output")
        return problems
    if summary["passed"] <= 0:
        problems.append(f"no tests passed (summary reports {summary['passed']} passed)")
    if summary["failed"]:
        problems.append(f"{summary['failed']} test(s) failed")
    if summary["error"]:
        problems.append(f"{summary['error']} collection/execution error(s)")
    return problems


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(f"usage: {Path(argv[0]).name} <pytest-output-file>", file=sys.stderr)
        return 2
    log_path = Path(argv[1])
    if not log_path.is_file():
        print(f"pytest output file not found: {log_path}", file=sys.stderr)
        return 2
    text = log_path.read_text(encoding="utf-8", errors="ignore")
    problems = check(text)
    summary = parse_summary(text)
    if summary is not None:
        print(
            "pytest summary: "
            f"passed={summary['passed']} failed={summary['failed']} "
            f"error={summary['error']} skipped={summary['skipped']}"
        )
    if problems:
        for problem in problems:
            print(f"FAILED: {problem}", file=sys.stderr)
        return 1
    print("pytest summary assertion: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
