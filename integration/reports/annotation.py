"""Postprocessing BLAST annotation of final backmt sequences, not humanization.

This module never changes sequences or reference databases. It reports the
best matching-chain hit in the configured human reconstructed-framework
database and explicitly distinguishes a closest match from a true assignment.
"""
from __future__ import annotations

import hashlib
import math
import re
import subprocess
import threading
from pathlib import Path

from .references import file_digest, catalog_entry, HUMAN_IMGT_FASTA_SHA256
from integration.paths import resolve_path

PARAMETERS = {"program": "blastp", "max_target_seqs": 50, "max_hsps": 1, "evalue": "1e-10"}
_CACHE = {}
_LOCK = threading.RLock()


def _identifier(value: str) -> str | None:
    match = re.fullmatch(r"(?:ref\|)?((?:VH|VK|VL)_\d+)\|?", value)
    return match[1] if match else None


def _snapshot(prefix: Path) -> dict:
    paths = {extension: Path(str(prefix) + extension) for extension in (".psq", ".pin", ".phr")}
    if not all(path.is_file() for path in paths.values()):
        raise ValueError("Configured human BLAST database is missing or incomplete.")
    snapshot = {extension: file_digest(path) for extension, path in paths.items()}
    for extension in (".pdb", ".pog", ".pos", ".pot", ".ptf", ".pto", ".pjs"):
        path = Path(str(prefix) + extension)
        if path.is_file():
            snapshot[extension] = file_digest(path)
    fasta = Path(str(prefix) + ".fasta")
    if fasta.is_file():
        snapshot[".fasta"] = file_digest(fasta)
    return snapshot


def _parse_hits(content: str, queries: dict, prefix: Path) -> dict:
    hits = {name: [] for name in queries}
    for line in content.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 13:
            raise ValueError("Malformed final-sequence BLAST output.")
        query, subject, identity, length, evalue, bits, qstart, qend, sstart, send, qseq, sseq, qlen = parts
        if query not in queries:
            raise ValueError("BLAST returned an unknown final-sequence query.")
        chain, sequence = queries[query]
        identifier = _identifier(subject)
        if not identifier or (chain == "VH" and not identifier.startswith("VH_")) or (
            chain == "VL" and not identifier.startswith(("VK_", "VL_"))
        ):
            continue
        identity, evalue, bits = float(identity), float(evalue), float(bits)
        length, qstart, qend, sstart, send, qlen = map(int, (length, qstart, qend, sstart, send, qlen))
        if (not all(math.isfinite(value) for value in (identity, evalue, bits))
                or not 0 <= identity <= 100 or evalue < 0 or bits < 0
                or length != len(qseq) or length != len(sseq) or qlen != len(sequence)
                or not (1 <= qstart <= qend <= len(sequence)) or not (1 <= sstart <= send)
                or qseq.replace("-", "") != sequence[qstart - 1:qend]):
            raise ValueError("Inconsistent BLAST annotation of a final sequence.")
        actual_identity = 100 * sum(left == right and left != "-" for left, right in zip(qseq, sseq)) / length
        if abs(identity - actual_identity) > 0.001:
            raise ValueError("Reported BLAST identity disagrees with the aligned residues.")
        hits[query].append({"identifier": identifier, "subject_id": subject, "identity": identity,
                            "length": length, "evalue": evalue, "bitscore": bits,
                            "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
                            "query_aligned": qseq, "subject_aligned": sseq})
    selected = {}
    for query, candidates in hits.items():
        if not candidates:
            selected[query] = {"Status": "No compatible-chain hit met the configured BLAST threshold."}
            continue
        # Keep BLAST's order when bit scores tie; disclose ambiguity rather than implying uniqueness.
        best = max(candidates, key=lambda item: item["bitscore"])
        entry = catalog_entry(prefix, best["identifier"])
        if (not entry or _identifier(entry["id"]) != best["identifier"]
                or entry["sequence"][best["sstart"] - 1:best["send"]] != best["subject_aligned"].replace("-", "")):
            raise ValueError("Selected final germline hit does not match the local reference entry.")
        best["description"] = entry["description"]
        best["genes"] = dict(re.findall(r"\b([VDJ]):([^\s|]+)", entry["description"]))
        best["equal_best_bitscore_hits_returned"] = sum(hit["bitscore"] == best["bitscore"] for hit in candidates)
        selected[query] = best
    return selected


def annotate_final_sequences(document: dict) -> dict:
    """Populate final annotations in the auxiliary document only; leave template data intact."""
    database = document["run"]["database_path"]
    if not database:
        for chain in document["chains"].values():
            for candidate in chain["candidates"]:
                candidate["germline_annotation"] = {"Status": "Final-sequence lookup could not run: human database path is not recorded."}
        document["final_germline_analysis"] = {"method": "blastp_on_final_backmt", "error": "Human database path is not recorded."}
        return document
    prefix = resolve_path(database)
    with _LOCK:
        try:
            snapshot = _snapshot(prefix)
            key_snapshot = (str(prefix), tuple(sorted(snapshot.items())))
            queries, missing, query_keys = {}, {}, {}
            for chain, data in document["chains"].items():
                for candidate in data["candidates"]:
                    sequence = candidate["backmutated_sequence"]
                    key = (key_snapshot, chain, sequence)
                    query_keys[(chain, candidate["number"])] = key
                    if key not in _CACHE and key not in missing.values():
                        query = f"final_{len(queries) + 1}"
                        queries[query] = (chain, sequence)
                        missing[query] = key
            if queries:
                fasta = "".join(f">{name}\n{sequence}\n" for name, (_, sequence) in queries.items())
                result = subprocess.run([
                    "blastp", "-db", str(prefix), "-outfmt",
                    "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send qseq sseq qlen",
                    "-max_target_seqs", str(PARAMETERS["max_target_seqs"]),
                    "-max_hsps", str(PARAMETERS["max_hsps"]), "-evalue", PARAMETERS["evalue"],
                ], input=fasta, text=True, capture_output=True, check=True, timeout=120)
                parsed = _parse_hits(result.stdout, queries, prefix)
                if snapshot != _snapshot(prefix):
                    raise ValueError("Human database changed while annotating final sequences.")
                for query, annotation in parsed.items():
                    _CACHE[missing[query]] = annotation
            origin = "IMGT-derived human reconstructed frameworks" if snapshot.get(".fasta") == HUMAN_IMGT_FASTA_SHA256 else "Human reference catalog; origin not verified"
            for chain, data in document["chains"].items():
                for candidate in data["candidates"]:
                    hit = _CACHE[query_keys[chain, candidate["number"]]]
                    if "Status" in hit:
                        candidate["germline_annotation"] = dict(hit)
                        continue
                    sequence = candidate["backmutated_sequence"]
                    candidate["final_germline_alignment"] = dict(hit)
                    candidate["germline_annotation"] = {
                        "Source": hit["subject_id"],
                        "Identity": f"{hit['identity']:.3f}% over {hit['length']} aligned positions",
                        "IMGT germline": " ".join(f"{kind}:{gene}" for kind, gene in hit["genes"].items()) or hit["description"],
                        "Query coverage": f"{100 * (hit['qend'] - hit['qstart'] + 1) / len(sequence):.2f}%",
                        "Bit score": str(hit["bitscore"]), "E-value": str(hit["evalue"]),
                        "Database origin": origin,
                        "Interpretation": "Closest matching-chain BLAST hit for the FINAL backmt sequence; not a definitive biological germline assignment and not necessarily the original grafting template.",
                        "Ties": f"{hit['equal_best_bitscore_hits_returned']} returned hits share the best bit score; additional ties may exist. First BLAST-ranked best-score hit selected.",
                        "Analysis scope": "Additional postprocessing against the database snapshot recorded here; humanization and Backmutation were not recalculated or changed.",
                    }
            document["final_germline_analysis"] = {
                "method": "blastp_on_final_backmt", "parameters": dict(PARAMETERS),
                "database_path": str(prefix), "database_sha256": snapshot, "database_origin": origin,
                "queries": {chain: {str(candidate["number"]): hashlib.sha256(candidate["backmutated_sequence"].encode()).hexdigest()
                                    for candidate in data["candidates"]} for chain, data in document["chains"].items()},
            }
        except (ValueError, OSError, subprocess.SubprocessError) as exc:
            for data in document["chains"].values():
                for candidate in data["candidates"]:
                    candidate["germline_annotation"] = {"Status": "Final-sequence germline lookup failed: " + str(exc)}
            document["final_germline_analysis"] = {"method": "blastp_on_final_backmt", "error": str(exc)}
    return document
