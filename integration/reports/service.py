"""Single entry point for result enrichment and organized TXT/JSON writing."""
from __future__ import annotations

import argparse
import json
import os
import tempfile
import threading
from pathlib import Path

from integration.paths import normalize_metadata, resolve_path, require_virtual_environment
from integration.fingerprints import add_comparison_marks
from .reader import JSON_NAME, TEXT_NAME, find_backmutation_report, build_organized_report
from .references import add_reference_info
from .annotation import annotate_final_sequences
from .rules import candidate_context, read_rule_source, load_rules_with_observation
from .render import render_text, compact_existing_text

_LOCK = threading.RLock()


def load_validated_rules(metadata: dict, reported: dict | None = None) -> dict:
    """Never treat an auxiliary organized TXT as an original rule log."""
    metadata = normalize_metadata(metadata)
    directory = Path(metadata["output_dir"])
    report = find_backmutation_report(directory)
    if report is None:
        return reported or {}
    lengths, _ = candidate_context(directory)
    original = read_rule_source(report, lengths)
    return load_rules_with_observation(metadata, reported=original)


def report_rule_data(document: dict) -> dict:
    results = {}
    for chain, data in document["chains"].items():
        for candidate in data["candidates"]:
            rules = [{**rule, "positions": [] if rule["positions"] == "-" else
                      [int(value) for value in rule["positions"].split(",")]}
                     for rule in candidate["rules"]]
            result = {"available": candidate["rules_available"], "source": candidate["rules_source"],
                      "rules": rules, "source_file": candidate.get("rules_source_file")}
            if candidate.get("rules_error"):
                result["error"] = candidate["rules_error"]
            results[(chain, candidate["number"])] = result
    return results


def build_report(directory: Path, metadata: dict | None = None) -> dict:
    directory = resolve_path(directory)
    metadata = normalize_metadata(metadata)
    raw = find_backmutation_report(directory)
    if raw is None:
        raise ValueError("No original Backmutation report.")
    document = build_organized_report(raw, directory / "humanized_vh.fasta", directory / "humanized_vl.fasta", metadata)
    document = add_reference_info(document, directory)
    if metadata:
        rules = load_validated_rules(metadata)
    else:
        lengths, _ = candidate_context(directory)
        rules = read_rule_source(raw, lengths)
    for chain, data in document["chains"].items():
        for candidate in data["candidates"]:
            evidence = rules.get((chain, candidate["number"]), {"available": False, "source": "missing_rules", "rules": []})
            candidate.update(rules_source=evidence["source"], rules_source_file=evidence.get("source_file"),
                             rules_available=evidence["available"],
                             rules=[{**rule, "positions": ",".join(map(str, rule["positions"])) or "-"}
                                    for rule in evidence["rules"]])
            if evidence.get("error"):
                candidate["rules_error"] = evidence["error"]
    document = annotate_final_sequences(document)
    document["complete_report_version"] = 1
    return add_comparison_marks(document, directory)


def _write_contents(contents: list[tuple[Path, str]]) -> None:
    if any(path.is_symlink() for path, _ in contents):
        raise ValueError("Refusing to overwrite a report symlink.")
    temporary = []
    try:
        for target, content in contents:
            fd, name = tempfile.mkstemp(prefix=".organized-report-", suffix=".tmp", dir=target.parent)
            temporary.append((Path(name), target))
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(content)
        for source, target in temporary:
            os.replace(source, target)
    finally:
        for source, _ in temporary:
            source.unlink(missing_ok=True)


def clean_saved_report(output_dir: Path) -> bool:
    """Normalize old saved TXT without executing any enrichment or scientific step."""
    with _LOCK:
        target = resolve_path(output_dir) / TEXT_NAME
        if target.is_symlink():
            raise ValueError("Refusing to overwrite a report symlink.")
        if not target.is_file():
            return False
        original = target.read_text(encoding="utf-8")
        compact = compact_existing_text(original)
        if compact == original:
            return False
        _write_contents([(target, compact)])
        return True


def ensure_report(output_dir: Path, metadata: dict | None = None) -> dict | None:
    with _LOCK:
        directory = resolve_path(output_dir)
        if find_backmutation_report(directory) is None:
            return None
        clean_saved_report(directory)
        document = build_report(directory, metadata)
        contents = [(directory / TEXT_NAME, render_text(document)),
                    (directory / JSON_NAME, json.dumps(document, indent=2, ensure_ascii=False) + "\n")]
        if not all(path.is_file() and not path.is_symlink()
                   and path.read_text(encoding="utf-8") == content for path, content in contents):
            _write_contents(contents)
        return document


def main() -> int:
    require_virtual_environment()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    try:
        directory = resolve_path(args.output_dir)
        path = directory.parent / "metadata.json"
        metadata = json.loads(path.read_text()) if path.is_file() else None
        if metadata and metadata.get("status") != "completed":
            raise ValueError("Execution has not completed.")
        result = ensure_report(directory, metadata)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Report unavailable: {exc}\n")
    if result:
        print(directory / TEXT_NAME)
        print(directory / JSON_NAME)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
