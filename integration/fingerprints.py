"""Stable comparison marks for results and recorded execution conditions."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

from integration.paths import HUMANIZE, BACK


def fingerprint(value: object) -> str:
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def record_execution(directory: Path, input_file: Path, database: Path,
                     mouse_database: Path, mode: str, level: str | None,
                     backmutation: bool) -> None:
    """Record current conditions at execution time, not retroactively for old runs."""
    def file_hash(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def database_hashes(prefix: Path) -> dict:
        return {extension: file_hash(Path(str(prefix) + extension))
                for extension in (".psq", ".pin", ".phr")
                if Path(str(prefix) + extension).is_file()}

    versions = {}
    for package in ("biopython", "anarcii", "torch"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    conditions = {
        "schema_version": 1,
        "input_sequences": input_file.read_text(encoding="utf-8").splitlines(),
        "mode": mode, "optimization_level": level, "backmutation": backmutation,
        "scientific_code": {path.name: file_hash(path) for path in
                            (HUMANIZE, HUMANIZE.with_name("cdr.py"),
                             HUMANIZE.with_name("optimizations.py"), BACK)},
        "human_database": database_hashes(database),
        "murine_database": database_hashes(mouse_database) if backmutation else None,
        "python": platform.python_version(), "packages": versions,
    }
    destination = directory / "execution_conditions.json"
    if destination.exists() or destination.is_symlink():
        raise ValueError("Execution conditions already exist; refusing to replace them.")
    with destination.open("x", encoding="utf-8") as handle:
        json.dump({"conditions": conditions, "sha256": fingerprint(conditions)}, handle,
                  indent=2, ensure_ascii=False)
        handle.write("\n")


def add_comparison_marks(document: dict, directory: Path) -> dict:
    result = {
        "schema_version": 1,
        "chains": {chain: [{key: candidate.get(key) for key in
                            ("number", "human_germline", "humanized_sequence",
                             "backmutated_sequence", "count", "changes", "rules_available", "rules")}
                           for candidate in data["candidates"]]
                   for chain, data in document["chains"].items()},
    }
    marks = {"scheme_version": 1, "result_sha256": fingerprint(result),
             "execution_sha256": None, "execution_status": "not_recorded"}
    saved = directory / "execution_conditions.json"
    if saved.is_symlink():
        raise ValueError("Execution conditions must not be a symlink.")
    if saved.is_file():
        execution = json.loads(saved.read_text(encoding="utf-8"))
        if execution["sha256"] != fingerprint(execution["conditions"]):
            raise ValueError("Execution condition fingerprint is inconsistent.")
        input_sequences = [document["chains"][chain]["input_sequence"] for chain in ("VH", "VL")]
        if execution["conditions"]["input_sequences"] != input_sequences:
            raise ValueError("Execution conditions refer to different input sequences.")
        marks.update(execution_sha256=execution["sha256"], execution_status="recorded_at_execution")
    document["comparison"] = marks
    return document
