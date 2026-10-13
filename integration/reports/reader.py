"""Read and validate original scientific outputs; never execute the pipeline."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

from Bio import SeqIO
from integration.paths import normalize_metadata

AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
RULES_UNAVAILABLE = "Rule information is not available for this original result."


def compare_sequences(before: str, after: str) -> list:
    if len(before) != len(after):
        raise ValueError(f"Different sequence lengths: {len(before)} and {len(after)}; comparison unavailable.")
    if not AA.fullmatch(before) or not AA.fullmatch(after):
        raise ValueError("Invalid or empty amino acid sequence; comparison unavailable.")
    return [
        {"position": i, "from": old, "to": new}
        for i, (old, new) in enumerate(zip(before, after), 1) if old != new
    ]


def load_humanized_candidates(output_dir: Path) -> dict:
    """Keep FASTA file order, regardless of header numbering."""
    from Bio import SeqIO

    result = {}
    for chain in ("VH", "VL"):
        path = Path(output_dir) / f"humanized_{chain.lower()}.fasta"
        if path.exists():
            for number, record in enumerate(SeqIO.parse(path, "fasta"), 1):
                sequence = str(record.seq)
                if not AA.fullmatch(sequence):
                    raise ValueError(f"Invalid FASTA sequence in {path.name}, candidate {number}.")
                result[(chain, number)] = sequence
    return result


def parse_original_backmutated_sequences(content: str) -> dict:
    result, counters, chain = {}, {}, None
    for line in content.splitlines():
        section = re.fullmatch(r"\s*Chain\s+(VH|VL)\s*", line)
        if section:
            chain = section[1]
            if chain in counters:
                raise ValueError(f"Duplicate Chain {chain} section; candidate order is ambiguous.")
            counters[chain] = 0
        elif re.match(r"\s*backmt\b", line):
            match = re.fullmatch(r"\s*backmt\s+([A-Z]+)\s*", line)
            if not chain or not match or not AA.fullmatch(match[1]):
                raise ValueError("Malformed original backmt line or missing chain section.")
            counters[chain] += 1
            result[(chain, counters[chain])] = match[1]
    if not result:
        raise ValueError("No backmt sequences found in original report.")
    return result


def parse_enriched_report(content: str) -> dict:
    result, key = {}, None
    for line in content.splitlines():
        text = line.strip()
        header = re.fullmatch(r"(VH|VL) candidate (\d+)", text)
        if header:
            key = (header[1], int(header[2]))
            if key[1] < 1 or key in result:
                raise ValueError("Duplicate or invalid candidate number in enriched report.")
            result[key] = {"count": None, "changes": [], "available": True, "format": "enriched"}
            continue
        if key is None:
            continue
        count = re.fullmatch(r"Back-mutations performed:\s*(\d+)", text)
        change = re.fullmatch(r"(\d+):\s*([A-Z])\s*->\s*([A-Z])", text)
        rule = re.fullmatch(r"rule\s+(.+?)\s+count=(\d+)\s+positions=(.*)", text)
        table_rule = re.fullmatch(r"(.+?)\s*\|\s*(\d+)\s*\|\s*(.+)", text)
        if count:
            if result[key]["count"] is not None:
                raise ValueError("Duplicate mutation count in enriched report.")
            result[key]["count"] = int(count[1])
        elif change:
            if int(change[1]) < 1 or change[2] == change[3] or not AA.fullmatch(change[2] + change[3]):
                raise ValueError("Invalid effective change in enriched report.")
            result[key]["changes"].append({"position": int(change[1]), "from": change[2], "to": change[3]})
        elif rule or table_rule:
            match = rule or table_rule
            result[key].setdefault("rules", []).append(
                {"name": match[1].strip(), "count": int(match[2]), "positions": match[3].strip()}
            )
        elif text.startswith("Rule |") or text == "Rules applied:":
            result[key].setdefault("rules", [])
        elif re.match(r"\d+:|Back-mutations performed:|rule\s", text):
            raise ValueError(f"Malformed enriched report line: {text}")
    if not result:
        raise ValueError("No candidates found in enriched report.")
    for key, data in result.items():
        positions = [change["position"] for change in data["changes"]]
        if data["count"] != len(positions) or len(set(positions)) != len(positions):
            raise ValueError(f"Inconsistent count or duplicate change positions for {key}.")
    return result


def find_backmutation_report(output_dir: Path) -> Path | None:
    directory = Path(output_dir)
    reports = [directory / "backmutation.txt", *sorted(directory.glob("*.backmutated"))]
    existing = [path for path in reports if path.is_file()]
    if len(existing) > 1:
        raise ValueError("Multiple backmutation reports found; candidate association is ambiguous.")
    return existing[0] if existing else None


def load_backmutation_results(output_dir: Path) -> dict:
    """Return explicit unavailable entries on format/association errors, never partial counts."""
    candidates = {}
    try:
        report = find_backmutation_report(output_dir)
        if report is None:
            return {}
        candidates = load_humanized_candidates(output_dir)
        content = report.read_text(encoding="utf-8")
        if re.search(r"^\s*(VH|VL) candidate \d+\s*$", content, re.M):
            result = parse_enriched_report(content)
            if not candidates or set(result) != set(candidates):
                raise ValueError("Report and FASTA candidate counts differ; association unavailable.")
            for key, data in result.items():
                for change in data["changes"]:
                    position = change["position"]
                    if position > len(candidates[key]) or candidates[key][position - 1] != change["from"]:
                        raise ValueError(f"Reported change disagrees with FASTA for {key} at {position}.")
        else:
            sequences = parse_original_backmutated_sequences(content)
            if not candidates or set(sequences) != set(candidates):
                raise ValueError("Report and FASTA candidate counts differ; association unavailable.")
            result = {}
            for key, sequence in sequences.items():
                changes = compare_sequences(candidates[key], sequence)
                result[key] = {
                    "count": len(changes), "changes": changes, "available": True,
                    "format": "original", "sequence": sequence,
                }
        for data in result.values():
            data["rules_available"] = "rules" in data
        return result
    except (ValueError, OSError) as exc:
        keys = candidates or {("VH", 1): "", ("VL", 1): ""}
        return {key: {"count": None, "changes": [], "available": False, "error": str(exc),
                      "rules_available": False} for key in keys}



SCHEMA_VERSION = 1
JSON_NAME = "backmutation_organized.json"
TEXT_NAME = "backmutation_organized.txt"
RULE_NAMES = {
    "murine_germline_difference": "Murine vs murine germline difference",
    "cysteine_proline": "Cysteine/Proline protection",
    "position_71": "Internal position 71 restoration",
    "fr4_motif": "FR4 motif correction",
}


def _candidate(number: int) -> dict:
    return {"number": number, "human_germline": None, "backmutated_sequence": None,
            "count": None, "changes": [], "rules": [], "rules_available": False,
            "germline_annotation": None}


def _rule(name: str, count: str, positions: str) -> dict:
    count = int(count)
    positions = positions.strip()
    if positions != "-" and not re.fullmatch(r"[1-9]\d*(?:\s*,\s*[1-9]\d*)*", positions):
        raise ValueError("Malformed rule positions: " + positions)
    numbers = [] if positions == "-" else [int(value) for value in positions.split(",")]
    if count != len(numbers) or len(set(numbers)) != len(numbers):
        raise ValueError("Rule count and positions disagree: " + name)
    return {"name": RULE_NAMES.get(name, name), "count": count,
            "positions": ",".join(map(str, numbers)) if numbers else "-"}


def parse_raw_report(content: str) -> dict:
    """Read original or organized text using explicit chain/candidate association."""
    chains, chain, active = {}, None, None
    enriched = bool(re.search(r"^\s*(VH|VL) candidate \d+\s*$", content, re.M))
    annotation = False
    for number, line in enumerate(content.splitlines(), 1):
        text = line.strip()
        section = re.fullmatch(r"CHAIN\s+(VH|VL)", text, re.I)
        header = re.fullmatch(r"(VH|VL) candidate ([1-9]\d*)", text)
        if section:
            name = section[1].upper()
            if name in chains:
                raise ValueError(f"Duplicate chain section {name} at line {number}.")
            chain = chains[name] = {"framework_map": None, "mouse_germline": None,
                                    "match": None, "input_sequence": None,
                                    "input_label": None, "candidates": []}
            active, annotation = None, False
            continue
        if header:
            name, index = header[1], int(header[2])
            if name not in chains:
                chains[name] = {"framework_map": None, "mouse_germline": None,
                                "match": None, "input_sequence": None,
                                "input_label": None, "candidates": []}
            chain = chains[name]
            if index != len(chain["candidates"]) + 1:
                raise ValueError(f"Out-of-order candidate at line {number}.")
            active = _candidate(index)
            chain["candidates"].append(active)
            annotation = False
            continue
        if not text:
            continue
        sequence_line = re.fullmatch(r"(hGerm|backmt|mGerm|input_)\s*:?[ \t]+([A-Z*-]+)", text)
        if sequence_line:
            label, sequence = sequence_line.groups()
            if chain is None or (label != "mGerm" and not AA.fullmatch(sequence)):
                raise ValueError(f"Invalid sequence or missing chain at line {number}.")
            if label == "hGerm":
                if not enriched:
                    active = _candidate(len(chain["candidates"]) + 1)
                    chain["candidates"].append(active)
                if active is None or active["human_germline"] is not None:
                    raise ValueError(f"Ambiguous hGerm candidate at line {number}.")
                active["human_germline"] = sequence
                annotation = False
            elif label == "backmt":
                if active is None or active["backmutated_sequence"] is not None:
                    raise ValueError(f"backmt without a unique active candidate at line {number}.")
                active["backmutated_sequence"] = sequence
            elif label == "mGerm":
                if chain["mouse_germline"] is not None:
                    raise ValueError("Duplicate mouse germline sequence.")
                chain["mouse_germline"] = sequence
            else:
                if chain["input_sequence"] is not None:
                    raise ValueError("Duplicate input sequence.")
                chain["input_sequence"], chain["input_label"] = sequence, label
            continue
        if re.match(r"^(hGerm|backmt|mGerm|input_)\b", text):
            raise ValueError(f"Malformed sequence line {number}: {text}")
        labelled = line.lstrip()
        match_prefix = next((prefix for prefix in ("match: ", "match  ", "match ")
                             if labelled.startswith(prefix)), None)
        if match_prefix and chain is not None:
            if chain["match"] is not None:
                raise ValueError("Duplicate mouse alignment.")
            chain["match"] = labelled[len(match_prefix):]
            continue
        if chain is not None and active is None and re.search(r"FR[1-4]-", line):
            chain["framework_map"] = line
            continue
        # Original back.py labels the input with the first six characters of its filename.
        named_input = re.fullmatch(r"(\S+)\s+([A-Z]+)", text)
        if chain is not None and active is None and named_input and AA.fullmatch(named_input[2]):
            if chain["input_sequence"] is not None:
                raise ValueError("Ambiguous input sequence label.")
            chain["input_label"] = named_input[1].rstrip(":")
            chain["input_sequence"] = named_input[2]
            continue
        rule = re.fullmatch(r"rule\s+(\S+)\s+count=(\d+)\s+positions=(\S+)", text)
        table_rule = re.fullmatch(r"(.+?)\s*\|\s*(\d+)\s*\|\s*(.+)", text)
        if rule or table_rule:
            if active is None or active["backmutated_sequence"] is None:
                raise ValueError("Rule without a completed active candidate.")
            fields = (rule or table_rule).groups()
            active["rules"].append(_rule(*fields))
            active["rules_available"] = True
            continue
        if text.startswith("rule "):
            raise ValueError(f"Malformed rule at line {number}.")
        if active is not None:
            count = re.fullmatch(r"Back-mutations performed:\s*(\d+)", text)
            change = re.fullmatch(r"([1-9]\d*):\s*([A-Z])\s*->\s*([A-Z])", text)
            if count:
                if active["count"] is not None:
                    raise ValueError("Duplicate reported mutation count.")
                active["count"] = int(count[1])
            elif change:
                active["changes"].append({"position": int(change[1]), "from": change[2], "to": change[3]})
            elif text.startswith("Back-mutations performed:") or re.match(r"\d+:", text):
                raise ValueError(f"Malformed mutation information at line {number}.")
            elif re.fullmatch(r"Rule\s*\|\s*Changes\s*\|\s*Positions", text):
                active["rules_available"] = True
            elif text == "Germline determined from this backmt sequence:":
                annotation = True
            elif annotation:
                field = re.fullmatch(r"(Source|Identity|IMGT germline):\s*(.+)", text)
                if field:
                    active["germline_annotation"] = active["germline_annotation"] or {}
                    active["germline_annotation"][field[1]] = field[2]
    if not chains or not any(chain["candidates"] for chain in chains.values()):
        raise ValueError("The Backmutation report is empty or contains no candidates.")
    for name, chain in chains.items():
        if not chain["candidates"]:
            raise ValueError(f"Chain {name} contains no complete candidates.")
        for candidate in chain["candidates"]:
            if candidate["backmutated_sequence"] is None:
                raise ValueError(f"Missing backmt for {name} candidate {candidate['number']}.")
            if not enriched and candidate["human_germline"] is None:
                raise ValueError(f"Missing hGerm for {name} candidate {candidate['number']}.")
    return {"chains": chains, "input_format": "organized" if enriched else "original"}


def _fingerprint(path: Path, directory: Path) -> dict:
    return {"path": str(path.relative_to(directory)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def build_organized_report(raw_report_path: Path, vh_fasta_path: Path,
                           vl_fasta_path: Path, metadata: dict | None = None) -> dict:
    """Pure read/validation step; no BLAST, sequence changes, or historical DB inference."""
    metadata = normalize_metadata(metadata)
    raw, vh, vl = [Path(p).resolve() for p in (raw_report_path, vh_fasta_path, vl_fasta_path)]
    if vh.parent != raw.parent or vl.parent != raw.parent:
        raise ValueError("Report and FASTA files must belong to the same output directory.")
    if metadata and metadata.get("output_dir") and Path(metadata["output_dir"]).resolve() != raw.parent:
        raise ValueError("Metadata and report belong to different executions.")
    document = parse_raw_report(raw.read_text(encoding="utf-8"))
    sources = [raw]
    for name, fasta in (("VH", vh), ("VL", vl)):
        records = list(SeqIO.parse(fasta, "fasta")) if fasta.is_file() else []
        if fasta.is_file():
            sources.append(fasta)
        chain = document["chains"].get(name)
        if len(records) != len(chain["candidates"] if chain else []):
            raise ValueError(f"Report/FASTA candidate counts differ for {name}.")
        if chain is None:
            continue
        for record, candidate in zip(records, chain["candidates"]):
            before, after = str(record.seq), candidate["backmutated_sequence"]
            changes = compare_sequences(before, after)
            if document["input_format"] == "organized":
                if candidate["count"] != len(changes) or candidate["changes"] != changes:
                    raise ValueError(f"Reported changes disagree with FASTA for {name} candidate {candidate['number']}.")
            candidate.update(humanized_sequence=before, proposal_id=record.id,
                             proposal_description=record.description, count=len(changes), changes=changes)
            for rule in candidate["rules"]:
                positions = [] if rule["positions"] == "-" else map(int, rule["positions"].split(","))
                if any(position > len(after) for position in positions):
                    raise ValueError("Reported rule position exceeds sequence length.")
        germ_file = raw.parent / ("hGerm" + name)
        if germ_file.is_file():
            sources.append(germ_file)
            germs = germ_file.read_text(encoding="utf-8").splitlines()
            if len(germs) != len(records):
                raise ValueError(f"hGerm/FASTA candidate counts differ for {name}.")
            for candidate, germ in zip(chain["candidates"], germs):
                if not AA.fullmatch(germ) or candidate["human_germline"] not in (None, germ):
                    raise ValueError(f"Human template differs between original artifacts for {name}.")
                candidate["human_germline"] = germ
        framework_file = raw.parent / ("FRfile" + name)
        if framework_file.is_file():
            sources.append(framework_file)
            chain["framework_sequences"] = framework_file.read_text(encoding="utf-8").splitlines()
    if metadata and metadata.get("input_file"):
        input_file = Path(metadata["input_file"]).resolve()
        if not input_file.is_relative_to(raw.parent.parent):
            raise ValueError("Input file does not belong to this execution directory.")
        if input_file.is_file():
            sequences = input_file.read_text(encoding="utf-8").splitlines()
            if len(sequences) != 2 or any(not AA.fullmatch(sequence) for sequence in sequences):
                raise ValueError("Invalid original VH/VL input file.")
            sources.append(input_file)
            for name, sequence in zip(("VH", "VL"), sequences):
                chain = document["chains"].get(name)
                if chain is not None:
                    if chain["input_sequence"] not in (None, sequence):
                        raise ValueError(f"Original input and report sequence differ for {name}.")
                    chain["input_sequence"] = sequence
    document["schema_version"] = SCHEMA_VERSION
    document["position_numbering"] = "one-based sequence indices (not Kabat/IMGT)"
    document["run"] = {key: (metadata or {}).get(key) for key in (
        "run_id", "work_name", "created_at", "run_mode", "optimization_level", "backmutation",
        "database_path", "murine_database_source", "murine_database_psq_sha256")}
    document["run"]["murine_database_source"] = document["run"]["murine_database_source"] or "Not recorded for this execution"
    document["sources"] = [_fingerprint(path, raw.parent.parent) for path in sources]
    document["totals"] = {
        "candidates": sum(len(chain["candidates"]) for chain in document["chains"].values()),
        "effective_changes": sum(candidate["count"] for chain in document["chains"].values()
                                 for candidate in chain["candidates"]),
    }
    return document



def dashboard_mutations(document: dict) -> dict:
    """Supply only the fields already shown by the existing dashboard."""
    results = {}
    for name, chain in document["chains"].items():
        for candidate in chain["candidates"]:
            data = {"count": candidate["count"], "changes": candidate["changes"],
                    "available": True, "rules_available": candidate["rules_available"],
                    "sequence": candidate["backmutated_sequence"]}
            if candidate["rules_available"]:
                data["rules"] = candidate["rules"]
            results[(name, candidate["number"])] = data
    return results
