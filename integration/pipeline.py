"""Orchestrate the unchanged original executables in an isolated run directory."""
from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path

from integration.reports.reader import load_backmutation_results, load_humanized_candidates

from integration.paths import PROJECT as PROJECT_DIR, PYTHON, HUMANIZE, BACK, MOUSE_DATABASE, resolve_path
from integration.fingerprints import record_execution


def build_humanize_command(input_file, output_dir, database_path, run_mode="complete", optimization_level=None):
    command = [str(PYTHON), str(HUMANIZE), str(Path(input_file).resolve()),
               "-o", str(Path(output_dir).resolve()), "--database", str(resolve_path(database_path))]
    if run_mode == "complete":
        level = str(optimization_level or "4")
        if level not in {"1", "2", "3", "4"}:
            raise ValueError("Optimization level must be 1, 2, 3 or 4.")
        command.extend(["--optimization", level])
    elif run_mode in {"graft", "graft_backmutation"}:
        command.append("--graft")
    elif run_mode in {"structures", "scores"}:
        command.append("--" + run_mode)
    else:
        raise ValueError(f"Invalid run mode: {run_mode}")
    return command


def execute_pipeline(input_file, output_dir, database_path, run_mode, optimization_level,
                     backmutation, log_handle) -> None:
    input_file, output_dir = Path(input_file).resolve(), Path(output_dir).resolve()
    do_back = backmutation or run_mode == "graft_backmutation"
    if do_back and run_mode not in {"complete", "graft", "graft_backmutation"}:
        raise ValueError("Backmutation requires complete or graft mode.")
    if not PYTHON.is_file():
        raise ValueError("Project .venv is missing.")
    output_dir.mkdir(parents=True, exist_ok=True)
    record_execution(output_dir, input_file, resolve_path(database_path), MOUSE_DATABASE,
                     run_mode, str(optimization_level or "4") if run_mode == "complete" else None, do_back)
    # humanize.py writes FRfile/hGerm in cwd. Never let concurrent runs share cwd.
    workspace = output_dir.parent / "pipeline_work"
    workspace.mkdir(exist_ok=False)
    env = dict(os.environ)
    env["PATH"] = str(PYTHON.parent) + os.pathsep + env.get("PATH", "")
    env["VIRTUAL_ENV"] = str(PYTHON.parent.parent)
    env["PYTHONUNBUFFERED"] = "1"

    def run(command, *, stdout=None, input_text=None):
        log_handle.write("COMMAND: " + shlex.join(command) + "\n")
        log_handle.flush()
        subprocess.run(command, cwd=workspace, env=env, input=input_text, text=True,
                       stdout=stdout if stdout is not None else log_handle,
                       stderr=log_handle, check=True)

    run(build_humanize_command(input_file, output_dir, database_path, run_mode, optimization_level))
    if run_mode in {"structures", "scores"}:
        return
    for name in ("FRfileVH", "FRfileVL", "hGermVH", "hGermVL"):
        source = workspace / name
        if source.exists():
            source.rename(output_dir / name)
    if not do_back:
        return

    candidates = load_humanized_candidates(output_dir)
    sequences = input_file.read_text(encoding="utf-8").splitlines()
    if len(sequences) != 2:
        raise ValueError("Backmutation input must have exactly two lines: VH and VL.")
    for chain, sequence in zip(("VH", "VL"), sequences):
        chain_candidates = [seq for (kind, _), seq in candidates.items() if kind == chain]
        if not chain_candidates:
            raise ValueError(f"No {chain} candidates; Backmutation cannot run.")
        # back.py indexes directly by input positions: refuse unsafe length mismatches.
        if any(len(candidate) != len(sequence) for candidate in chain_candidates):
            raise ValueError(f"{chain} candidate/input lengths differ; original Backmutation unavailable.")
        for prefix in ("FRfile", "hGerm"):
            if not (output_dir / (prefix + chain)).is_file():
                raise ValueError(f"Missing {prefix + chain}; Backmutation cannot run.")
        germs = (output_dir / ("hGerm" + chain)).read_text().splitlines()
        if len(germs) != len(chain_candidates):
            raise ValueError(f"{chain} hGerm/FASTA candidate counts differ; association unavailable.")
        xml = output_dir / (chain + ".xml")
        with xml.open("w", encoding="utf-8") as handle:
            run(["blastp", "-db", str(MOUSE_DATABASE), "-max_target_seqs", "1",
                 "-max_hsps", "1", "-outfmt", "5"], stdout=handle, input_text=sequence + "\n")
        from Bio.Blast import NCBIXML
        with xml.open() as handle:
            record = NCBIXML.read(handle)
        if not record.alignments or not record.alignments[0].hsps or not record.alignments[0].hsps[0].query:
            raise ValueError(f"No valid murine BLAST alignment for {chain}.")
    # Report is raw stdout from the original back.py; never enrich/replace it.
    report = output_dir / (input_file.stem + ".backmutated")
    with report.open("w", encoding="utf-8") as handle:
        run([str(PYTHON), str(BACK), str(input_file), str(output_dir),
             str(output_dir / "VH.xml"), str(output_dir / "VL.xml")], stdout=handle)
    results = load_backmutation_results(output_dir)
    if not results or any(not result["available"] for result in results.values()):
        raise ValueError("Backmutation report is incomplete or inconsistent; see the auxiliary report.")
