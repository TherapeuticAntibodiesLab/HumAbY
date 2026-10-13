"""Run the unchanged scientific pipeline from the organized project layout."""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import traceback
from datetime import datetime
from pathlib import Path

from integration.paths import PROJECT, HUMAN_DATABASE, require_virtual_environment, resolve_path
from integration.pipeline import execute_pipeline
from integration.reports.service import ensure_report


def main() -> int:
    require_virtual_environment()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Two-line VH/VL input file")
    parser.add_argument("--mode", choices=("complete", "graft", "graft_backmutation", "structures", "scores"), default="complete")
    parser.add_argument("--optimization", choices=("1", "2", "3", "4"), default="4")
    parser.add_argument("--backmutation", action="store_true")
    parser.add_argument("--database", type=Path, default=HUMAN_DATABASE)
    parser.add_argument("--output", type=Path, help="New run directory; must not already exist")
    args = parser.parse_args()
    original = resolve_path(args.input)
    if not original.is_file():
        parser.error(f"Input file not found: {original}")
    backmutation = args.backmutation or args.mode == "graft_backmutation"
    if backmutation and args.mode not in {"complete", "graft", "graft_backmutation"}:
        parser.error("Backmutation requires complete or graft mode.")
    now = datetime.now()
    name = original.stem + "_" + now.strftime("%Y%m%d_%H%M%S_%f")
    directory = resolve_path(args.output) if args.output else PROJECT / "results" / "manual" / name
    if directory.exists() or directory.is_symlink():
        parser.error(f"Output directory already exists; it will not be overwritten: {directory}")
    directory.mkdir(parents=True)
    for part in ("input", "output", "logs"):
        (directory / part).mkdir()
    input_file = directory / "input" / original.name
    shutil.copyfile(original, input_file)
    output = directory / "output"
    log = directory / "logs" / "pipeline.log"
    metadata = {
        "run_id": directory.name, "work_name": original.stem,
        "created_at": now.isoformat(timespec="seconds"), "status": "running",
        "input_file": str(input_file), "output_dir": str(output), "log_file": str(log),
        "database_path": str(resolve_path(args.database)), "run_mode": args.mode,
        "optimization_level": args.optimization if args.mode == "complete" else None,
        "backmutation": backmutation,
    }
    metadata_path = directory / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    try:
        with log.open("w", encoding="utf-8") as handle:
            execute_pipeline(input_file, output, args.database, args.mode,
                             metadata["optimization_level"], backmutation, handle)
        metadata.update(status="completed", return_code=0)
    except Exception as exc:
        with log.open("a", encoding="utf-8") as handle:
            handle.write("\n" + traceback.format_exc())
        metadata.update(status="failed", return_code=1, error=str(exc))
    metadata["finished_at"] = datetime.now().isoformat(timespec="seconds")
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if metadata["status"] == "failed":
        print(f"Pipeline failed. Details: {log}", file=sys.stderr)
        return 1
    try:
        ensure_report(output, metadata)
    except Exception as exc:
        print(f"Scientific results saved; auxiliary report unavailable: {exc}", file=sys.stderr)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
