"""Validate reported rules or observe an exact-output replay of unchanged back.py.

Rule parsing is read-only; verified observations are written only to auxiliary
sidecars. The original scientific rules are never copied or reimplemented.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from Bio import SeqIO

import os
import subprocess
import tempfile
import threading
from datetime import datetime

from integration.paths import PROJECT, PYTHON, normalize_metadata
from .reader import find_backmutation_report, parse_original_backmutated_sequences, RULE_NAMES as RULE_LABELS
from .observer_worker import CORE, RULE_IDS, assignment_sites

WORKER = Path(__file__).with_name("observer_worker.py")

_IDS = {label: identifier for identifier, label in RULE_LABELS.items()}
RAW_RULE_PATTERN = re.compile(r"^rule\s+(\S+)\s+count=(\d+)\s+positions=(\S+)$")
TABLE_RULE_PATTERN = re.compile(r"^(.+?)\s*\|\s*(\d+)\s*\|\s*(.+)$")
TABLE_HEADER_PATTERN = re.compile(r"^Rule\s*\|\s*Changes\s*\|\s*Positions$")
CANDIDATE_PATTERN = re.compile(r"^(VH|VL) candidate ([1-9]\d*)$")
SCHEMA_VERSION = 1


def _missing() -> dict:
    return {"available": False, "source": "missing_rules", "rules": []}


def _parse_rule(name: str, count_text: str, positions_text: str) -> dict:
    count = int(count_text)
    if positions_text == "-":
        positions = []
    elif re.fullmatch(r"[1-9]\d*(?:\s*,\s*[1-9]\d*)*", positions_text):
        positions = [int(part) for part in positions_text.split(",")]
    else:
        raise ValueError(f"Invalid positions for rule {name}: {positions_text}")
    if count != len(set(positions)) or len(set(positions)) != len(positions):
        raise ValueError(f"Rule count does not match distinct positions: {name}")
    identifier = _IDS.get(name.strip(), name.strip())
    return {"id": identifier, "name": RULE_LABELS.get(identifier, identifier),
            "count": count, "positions": positions}


def _parse(content: str) -> dict:
    result, counters, sections = {}, {"VH": 0, "VL": 0}, set()
    chain, key = None, None
    in_rule_section = False
    explicit_candidates = bool(re.search(r"^\s*(VH|VL) candidate \d+\s*$", content, re.M))
    for number, line in enumerate(content.splitlines(), 1):
        text = line.strip()
        section = re.fullmatch(r"CHAIN\s+(VH|VL)", text, re.I)
        header = CANDIDATE_PATTERN.fullmatch(text)
        if section:
            chain = section[1].upper()
            if chain in sections:
                raise ValueError(f"Duplicate {chain} section at line {number}.")
            sections.add(chain)
            key = None  # Never carry a VH candidate into the VL section.
            in_rule_section = False
            continue
        if header:
            chain, index = header[1], int(header[2])
            if index != counters[chain] + 1:
                raise ValueError(f"Duplicate or out-of-order candidate at line {number}.")
            counters[chain] = index
            key = (chain, index)
            result[key] = _missing()
            in_rule_section = False
            continue
        if re.match(r"^hGerm(?:\s|:)", text) and not explicit_candidates:
            if chain is None:
                raise ValueError("hGerm appears before a chain section.")
            counters[chain] += 1
            key = (chain, counters[chain])
            result[key] = _missing()
            in_rule_section = False
            continue
        if text == "Rules applied:":
            if key is None:
                raise ValueError("Rules section without an active candidate.")
            in_rule_section = True
            continue
        table_header = TABLE_HEADER_PATTERN.fullmatch(text)
        raw = RAW_RULE_PATTERN.fullmatch(text)
        row = TABLE_RULE_PATTERN.fullmatch(text)
        if table_header or raw or row:
            if key is None:
                raise ValueError(f"Rule information before an active candidate at line {number}.")
            entry = result[key]
            entry.update(available=True, source="reported_empty")
            if raw or row:
                rule = _parse_rule(*(raw or row).groups())
                if any(previous["id"] == rule["id"] for previous in entry["rules"]):
                    raise ValueError("Duplicate rule type for one candidate: " + rule["id"])
                entry["rules"].append(rule)
            continue
        if re.match(r"^rule\b", text) or ("|" in text and key is not None
                                         and (in_rule_section or result[key]["available"])):
            raise ValueError(f"Malformed rule information at line {number}: {text}")
        if text == "No rule produced residue changes.":
            if key is None:
                raise ValueError("Zero-rule declaration without a candidate.")
            result[key].update(available=True, source="reported_empty")
        elif "Rule information is not available" in text:
            if key is not None and result[key]["available"]:
                raise ValueError("Report contradicts itself about rule availability.")
    for entry in result.values():
        if entry["available"]:
            entry["source"] = "reported" if any(rule["count"] > 0 for rule in entry["rules"]) else "reported_empty"
    return result


def parse_raw_rules(content: str) -> dict:
    return _parse(content)


def parse_organized_rules(content: str) -> dict:
    return _parse(content)


def candidate_context(output_dir: Path) -> tuple[dict, dict]:
    lengths, hashes = {}, {}
    for chain in ("VH", "VL"):
        path = output_dir / f"humanized_{chain.lower()}.fasta"
        if path.is_file():
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            for number, record in enumerate(SeqIO.parse(path, "fasta"), 1):
                lengths[(chain, number)] = len(record.seq)
    return lengths, hashes


def read_rule_source(path: Path, lengths: dict) -> dict:
    result = _parse(path.read_text(encoding="utf-8"))
    if lengths and result and set(result) != set(lengths):
        raise ValueError("Rule report candidates do not match the VH/VL FASTA order and counts.")
    for key, entry in result.items():
        for rule in entry["rules"]:
            if key in lengths and any(position > lengths[key] for position in rule["positions"]):
                raise ValueError(f"Rule position exceeds the sequence length for {key}.")
        entry["source_file"] = path.name
    return result


def _invalid(keys, error: Exception, source: str | None = None) -> dict:
    return {key: {"available": False, "source": "invalid", "rules": [],
                  "error": str(error), "source_file": source} for key in keys}


def load_candidate_rules(output_dir: Path, explicit_sources=()) -> dict:
    """Select one complete source per run; ignore stale/unverified JSON caches."""
    directory = Path(output_dir).resolve()
    try:
        lengths, fasta_hashes = candidate_context(directory)
    except (ValueError, OSError) as exc:
        return _invalid([("VH", 1), ("VL", 1)], exc)
    keys = list(lengths)
    cached_file = directory / "backmutation_rules.json"
    try:
        if cached_file.is_file():
            cached = json.loads(cached_file.read_text(encoding="utf-8"))
            source = (directory.parent / cached["source_file"]).resolve()
            if (cached.get("schema_version") == SCHEMA_VERSION
                    and source.is_relative_to(directory.parent) and source.is_file()
                    and cached.get("source_sha256") == hashlib.sha256(source.read_bytes()).hexdigest()
                    and cached.get("fasta_sha256") == fasta_hashes):
                parsed = read_rule_source(source, lengths)
                rows = [{"chain": key[0], "candidate": key[1], **entry} for key, entry in parsed.items()]
                # A matching hash alone cannot validate fabricated/inferred cache rules.
                if cached.get("candidates") == rows:
                    return parsed
    except (OSError, ValueError, KeyError, TypeError):
        pass

    sources = [directory / "backmutation.txt", *sorted(directory.glob("*.backmutated")),
               directory / "backmutation_organized.txt"]
    for path in explicit_sources:
        path = Path(path).resolve()
        if not path.is_relative_to(directory.parent):
            return _invalid(keys, ValueError("Explicit rule source belongs to another execution."))
        sources.append(path)
    missing = {key: _missing() for key in keys}
    seen = set()
    for source in sources:
        if not source.is_file() or source.resolve() in seen:
            continue
        if not source.resolve().is_relative_to(directory.parent):
            return _invalid(keys, ValueError("Rule source resolves outside this execution."), source.name)
        seen.add(source.resolve())
        try:
            parsed = read_rule_source(source, lengths)
        except (ValueError, OSError) as exc:
            return _invalid(keys or [("VH", 1), ("VL", 1)], exc, source.name)
        if parsed:
            missing = parsed
        if any(entry["available"] for entry in parsed.values()):
            return parsed  # Do not fill remaining candidates from a different source.
    return missing


def write_rules_sidecar(data: dict, destination: Path, source_file: Path) -> None:
    """Optional derived cache with source/FASTA hashes; never overwrites a source."""
    if Path(destination).is_symlink():
        raise ValueError("Refusing to overwrite a sidecar symlink.")
    destination, source_file = Path(destination).resolve(), Path(source_file).resolve()
    directory = destination.parent
    if destination.name != "backmutation_rules.json" or destination == source_file:
        raise ValueError("Only an additional backmutation_rules.json may be written.")
    if not source_file.is_relative_to(directory.parent):
        raise ValueError("Rule source belongs to another execution.")
    lengths, fasta_hashes = candidate_context(directory)
    parsed = read_rule_source(source_file, lengths)
    if data != parsed:
        raise ValueError("Sidecar rules must exactly match explicitly reported source data.")
    payload = {"schema_version": SCHEMA_VERSION,
               "source_file": str(source_file.relative_to(directory.parent)),
               "source_sha256": hashlib.sha256(source_file.read_bytes()).hexdigest(),
               "fasta_sha256": fasta_hashes,
               "candidates": [{"chain": key[0], "candidate": key[1], **entry} for key, entry in parsed.items()]}
    original_target = directory / "backmutation_rules.json"
    if original_target.is_symlink():
        raise ValueError("Refusing to overwrite a sidecar symlink.")
    fd, temporary = tempfile.mkstemp(prefix=".backmutation-rules-", suffix=".tmp", dir=directory)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        os.replace(temporary, destination)
    finally:
        Path(temporary).unlink(missing_ok=True)



def merge_rules(mutations: dict, rules_by_candidate: dict) -> dict:
    """Add presentation fields while preserving all scientific counts/changes."""
    merged = {}
    for key, data in mutations.items():
        value = dict(data)
        reported = rules_by_candidate.get(key, _missing())
        value.pop("rules", None)
        value.pop("rules_error", None)
        value.update(rules_available=reported["available"], rules_source=reported["source"],
                     rules_source_file=reported.get("source_file"))
        if reported["source"] == "invalid":
            value["rules_error"] = reported["error"]
        elif reported["available"]:
            value["rules"] = [{**rule, "positions": ",".join(map(str, rule["positions"])) or "-"}
                              for rule in reported["rules"]]
        merged[key] = value
    return merged


JSON_NAME = "backmutation_rules_observed.json"
TEXT_NAME = "backmutation_rules_observed.txt"
_LOCK = threading.RLock()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inputs(metadata: dict) -> tuple[Path, Path, Path, dict]:
    metadata = normalize_metadata(metadata)
    directory = Path(metadata["output_dir"]).resolve()
    input_file = Path(metadata["input_file"]).resolve()
    if not input_file.is_relative_to(directory.parent):
        raise ValueError("Original input belongs to another execution.")
    report = find_backmutation_report(directory)
    if report is None:
        raise ValueError("No original Backmutation report.")
    files = [input_file, report]
    files += [directory / name for name in ("humanized_vh.fasta", "humanized_vl.fasta",
              "FRfileVH", "FRfileVL", "hGermVH", "hGermVL", "VH.xml", "VL.xml")]
    if any(not path.is_file() for path in files):
        raise ValueError("Frozen input/FASTA/framework/germline/XML artifacts are incomplete.")
    if any(not path.resolve().is_relative_to(directory.parent) for path in files):
        raise ValueError("An observation artifact resolves to another execution.")
    sources = {str(path.relative_to(directory.parent)): _sha(path) for path in files}
    return directory, input_file, report, sources


def _validate_candidates(payload: dict, directory: Path, report: Path) -> dict:
    sites, _ = assignment_sites(CORE)
    if payload["core_sha256"] != _sha(CORE):
        raise ValueError("Observed core hash is different.")
    observed_sites = {int(line): rule for line, rule in payload["sites"].items()}
    if observed_sites != sites:
        raise ValueError("Observation refers to different assignment sites.")
    proposals = {}
    for chain in ("VH", "VL"):
        for index, record in enumerate(SeqIO.parse(directory / f"humanized_{chain.lower()}.fasta", "fasta"), 1):
            proposals[(chain, index)] = (record.id, str(record.seq))
    final = parse_original_backmutated_sequences(report.read_text(encoding="utf-8"))
    observed = {}
    for candidate in payload["candidates"]:
        key = (candidate["chain"], candidate["candidate"])
        if key in observed or key not in proposals or key not in final:
            raise ValueError("Observed candidate association is invalid.")
        identifier, sequence = proposals[key]
        if candidate["proposal_id"] != identifier or candidate["before_sequence"] != sequence:
            raise ValueError("Observed candidate differs from its original FASTA entry.")
        if candidate["after_sequence"] != final[key]:
            raise ValueError("Observed output differs from the saved scientific result.")
        state = list(sequence)
        positions = {rule: set() for rule in RULE_IDS}
        for event in candidate["events"]:
            rule, line, position = event["rule_id"], event["core_line"], event["position"]
            if sites.get(line) != rule or rule not in positions or not (1 <= position <= len(state)):
                raise ValueError("Observed event has an invalid rule/site/position.")
            if state[position - 1] != event["from"] or event["from"] == event["to"]:
                raise ValueError("Observed residue transition is inconsistent.")
            state[position - 1] = event["to"]
            positions[rule].add(position)
        if "".join(state) != final[key]:
            raise ValueError("Recorded assignment events do not reproduce the saved final sequence.")
        rules = [{"id": rule, "name": RULE_LABELS[rule], "count": len(positions[rule]),
                  "positions": sorted(positions[rule])} for rule in RULE_IDS]
        observed[key] = {"available": True, "source": "observed_replay", "rules": rules,
                         "source_file": JSON_NAME}
    if set(observed) != set(proposals) or set(observed) != set(final):
        raise ValueError("Observation is incomplete.")
    return observed


def _render(payload: dict, rules: dict) -> str:
    lines = ["BACKMUTATION RULES — VERIFIED REPLAY", "===================================", "",
             "The original back.py was executed without changing its functions or duplicating its algorithm.",
             "This is an external observation of a replay, not a rule log emitted by the historical run.",
             "Replay stdout and final sequences exactly match the saved original report.",
             "Positions are one-based sequence indices, not Kabat/IMGT numbering.",
             "Rule counts measure changed positions during actual assignment execution; counts can overlap.",
             "Core SHA-256: " + payload["core_sha256"], ""]
    for key, data in rules.items():
        lines += [f"{key[0]} candidate {key[1]}", "Rules recorded by external execution observation:",
                  "  Rule | Changes | Positions"]
        lines += [f"  {rule['name']} | {rule['count']} | {','.join(map(str, rule['positions'])) or '-'}"
                  for rule in data["rules"]]
        lines.append("")
    lines += ["FROZEN SOURCE FILES (SHA-256)"]
    lines += [f"{path}: {digest}" for path, digest in payload["source_sha256"].items()]
    return "\n".join(lines) + "\n"


def _write(directory: Path, payload: dict, rules: dict) -> None:
    contents = [(directory / TEXT_NAME, _render(payload, rules)),
                (directory / JSON_NAME, json.dumps(payload, indent=2, ensure_ascii=False) + "\n")]
    temporary = []
    try:
        for target, text in contents:
            if target.is_symlink():
                raise ValueError("Refusing to overwrite an observation symlink.")
            fd, name = tempfile.mkstemp(prefix=".rule-observation-", suffix=".tmp", dir=directory)
            temporary.append((Path(name), target))
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(text)
        for source, target in temporary:
            os.replace(source, target)
    finally:
        for source, _ in temporary:
            source.unlink(missing_ok=True)


def ensure_observation(metadata: dict) -> dict:
    with _LOCK:
        directory, input_file, report, sources = _inputs(metadata)
        cache = directory / JSON_NAME
        try:
            payload = json.loads(cache.read_text(encoding="utf-8"))
            if (isinstance(payload, dict) and payload.get("schema_version") == 1 and payload.get("source") == "observed_replay"
                    and payload.get("source_sha256") == sources
                    and payload.get("worker_sha256") == _sha(WORKER)
                    and payload.get("controller_sha256") == _sha(Path(__file__))
                    and payload.get("stdout_sha256") == _sha(report)):
                rules = _validate_candidates(payload, directory, report)
                if (directory / TEXT_NAME).is_file() and (directory / TEXT_NAME).read_text() == _render(payload, rules):
                    return rules
                _write(directory, payload, rules)
                return rules
        except (OSError, ValueError, KeyError, TypeError):
            pass
        if not PYTHON.is_file():
            raise ValueError("Project virtual environment is missing.")
        with tempfile.TemporaryDirectory(prefix="backmutation-observation-") as temporary:
            observation = Path(temporary) / "observed.json"
            process = subprocess.run([str(PYTHON), "-B", "-m", "integration.reports.observer_worker", str(input_file), str(directory),
                                      str(directory / "VH.xml"), str(directory / "VL.xml"), str(observation)],
                                     capture_output=True, check=True, timeout=60,
                                     cwd=str(PROJECT))
            if process.stdout != report.read_bytes():
                raise ValueError("Observed replay stdout differs from the saved original report; rules were not published.")
            payload = json.loads(observation.read_text(encoding="utf-8"))
        _, _, _, after_sources = _inputs(metadata)
        if sources != after_sources:
            raise ValueError("Original artifacts changed during observation.")
        payload.update(schema_version=1, source="observed_replay", verification="exact_stdout_and_final_sequences",
                       recorded_at=datetime.now().isoformat(timespec="seconds"),
                       worker_sha256=_sha(WORKER), stdout_sha256=_sha(report), source_sha256=sources,
                       controller_sha256=_sha(Path(__file__)),
                       historical_claim="This records the current unchanged core replay, not an original historical rule log.")
        rules = _validate_candidates(payload, directory, report)
        _write(directory, payload, rules)
        return rules


def load_rules_with_observation(metadata: dict, reported: dict | None = None) -> dict:
    """Prefer original explicit logs; observe complete original outputs on demand."""
    metadata = normalize_metadata(metadata)
    directory = Path(metadata["output_dir"])
    reported = load_candidate_rules(directory) if reported is None else reported
    if any(data.get("available") or data.get("source") == "invalid" for data in reported.values()):
        return reported
    if metadata.get("status") != "completed" or not reported:
        return reported
    if "input_file" not in metadata or not Path(metadata["input_file"]).is_file():
        return reported
    try:
        return ensure_observation(metadata)
    except (ValueError, OSError, subprocess.SubprocessError, KeyError) as exc:
        return {key: {"available": False, "source": "invalid", "rules": [],
                      "error": "Verified rule observation unavailable: " + str(exc)} for key in reported}
