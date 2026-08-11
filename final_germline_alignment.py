#!/usr/bin/env python3
"""Post-process final backmutated sequences against the human germline DB.

This module is deliberately independent from the humanization pipeline.  Its
public entry point receives the raw output produced by ``back/back.py`` and
creates its own audit log and human-readable report.
"""

from __future__ import annotations

import logging
import re
import subprocess
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import DefaultDict, Dict, List, Sequence, Tuple


VALID_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
REPORT_NAME = "final_germline_top5.txt"
LOG_NAME = "final_germline_alignment.log"
BLAST_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
    "evalue bitscore qseq sseq stitle"
)


class FinalGermlineAlignmentError(RuntimeError):
    """Raised when final-sequence germline analysis cannot be completed."""


def extract_final_sequences(raw_output: str) -> List[Tuple[str, str]]:
    """Extract ordered ``(query_id, sequence)`` pairs from back.py output."""
    chain = None
    counts: DefaultDict[str, int] = defaultdict(int)
    sequences: List[Tuple[str, str]] = []

    for raw_line in raw_output.splitlines():
        line = raw_line.strip()
        chain_match = re.match(r"^Chain\s+(VH|VL)\b", line, re.IGNORECASE)
        if chain_match:
            chain = chain_match.group(1).upper()
            continue
        if not line.startswith("backmt"):
            continue
        if chain is None:
            raise FinalGermlineAlignmentError(
                "A backmt sequence was found before its VH/VL chain header"
            )
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            raise FinalGermlineAlignmentError("Empty backmt sequence in backmutation output")
        sequence = parts[1].strip().upper()
        invalid = sorted(set(sequence) - VALID_AMINO_ACIDS)
        if not sequence or invalid:
            raise FinalGermlineAlignmentError(
                f"Invalid {chain} backmt sequence; unsupported residues={invalid}"
            )
        counts[chain] += 1
        sequences.append((f"{chain}_final_{counts[chain]}", sequence))

    if not sequences:
        raise FinalGermlineAlignmentError("No final backmt sequences found")
    return sequences


def _midline(query: str, subject: str) -> str:
    return "".join("|" if q == s else " " for q, s in zip(query, subject))


def _parse_hits(stdout: str) -> Dict[str, List[Dict[str, object]]]:
    hits: DefaultDict[str, List[Dict[str, object]]] = defaultdict(list)
    for line_number, line in enumerate(stdout.splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("\t", 14)
        if len(fields) != 15:
            raise FinalGermlineAlignmentError(
                f"Unexpected BLAST row at line {line_number}: expected 15 fields"
            )
        (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,
         sstart, send, evalue, bitscore, qseq, sseq, title) = fields
        hits[qseqid].append({
            "sseqid": sseqid,
            "pident": float(pident),
            "length": int(length),
            "mismatch": int(mismatch),
            "gapopen": int(gapopen),
            "qstart": int(qstart),
            "qend": int(qend),
            "sstart": int(sstart),
            "send": int(send),
            "evalue": float(evalue),
            "bitscore": float(bitscore),
            "qseq": qseq,
            "sseq": sseq,
            "title": title,
        })
    return dict(hits)


def _write_report(
    path: Path,
    database: Path,
    sequences: Sequence[Tuple[str, str]],
    hits: Dict[str, List[Dict[str, object]]],
    top_n: int,
) -> None:
    lines = [
        "FINAL BACKMUTATED SEQUENCES VS HUMAN GERMLINES",
        "=" * 55,
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        f"Human germline database: {database}",
        f"Ranking: BLAST bit score descending (top {top_n} per final sequence)",
        "Identity and alignment coordinates refer to the local BLAST alignment.",
        "",
    ]
    for query_id, sequence in sequences:
        query_hits = sorted(
            hits.get(query_id, []),
            key=lambda hit: (-float(hit["bitscore"]), float(hit["evalue"])),
        )[:top_n]
        lines.extend([
            f"QUERY: {query_id}",
            "-" * 55,
            f"Final sequence ({len(sequence)} aa): {sequence}",
            f"Human germline matches returned: {len(query_hits)}",
            "",
        ])
        if not query_hits:
            lines.extend(["No human germline match found.", ""])
            continue
        for rank, hit in enumerate(query_hits, 1):
            lines.extend([
                f"#{rank}  {hit['sseqid']}",
                f"Description: {hit['title'] or 'not available'}",
                (
                    f"Identity: {hit['pident']:.2f}% | Bit score: {hit['bitscore']:.1f} | "
                    f"E-value: {hit['evalue']:.3g} | Aligned length: {hit['length']} aa | "
                    f"Mismatches: {hit['mismatch']} | Gap openings: {hit['gapopen']}"
                ),
                f"Coordinates: query {hit['qstart']}-{hit['qend']}; germline {hit['sstart']}-{hit['send']}",
                f"Query     {hit['qseq']}",
                f"          {_midline(str(hit['qseq']), str(hit['sseq']))}",
                f"Germline  {hit['sseq']}",
                "",
            ])
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def run_final_germline_alignment(
    raw_backmutation_output: str,
    human_database: str,
    output_dir: str,
    top_n: int = 5,
) -> Path:
    """Align every final backmutated sequence and write a top-N report and log."""
    if top_n < 1:
        raise ValueError("top_n must be at least 1")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    report_path = output_path / REPORT_NAME
    log_path = output_path / LOG_NAME
    database = Path(human_database)

    audit = logging.getLogger(f"final_germline_alignment.{id(report_path)}")
    audit.setLevel(logging.INFO)
    audit.propagate = False
    handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    audit.addHandler(handler)

    query_path = None
    try:
        audit.info("FINAL GERMLINE ALIGNMENT START")
        audit.info("database=%s output_dir=%s top_n=%d", database, output_path, top_n)
        sequences = extract_final_sequences(raw_backmutation_output)
        audit.info("final sequences extracted=%d ids=%s", len(sequences), [x[0] for x in sequences])

        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False, encoding="utf-8") as query:
            query_path = Path(query.name)
            for query_id, sequence in sequences:
                query.write(f">{query_id}\n{sequence}\n")

        cmd = [
            "blastp", "-query", str(query_path), "-db", str(database),
            "-outfmt", f"6 {BLAST_FIELDS}", "-max_target_seqs", str(top_n),
            "-evalue", "1e-10",
        ]
        audit.info("running command=%s", " ".join(cmd))
        completed = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=120)
        if completed.stderr.strip():
            audit.warning("blastp stderr=%s", completed.stderr.strip())
        hits = _parse_hits(completed.stdout)
        for query_id, _ in sequences:
            audit.info("query=%s hits=%d", query_id, len(hits.get(query_id, [])))
        _write_report(report_path, database, sequences, hits, top_n)
        audit.info("report written=%s bytes=%d", report_path, report_path.stat().st_size)
        audit.info("FINAL GERMLINE ALIGNMENT END status=success")
        return report_path
    except Exception:
        audit.exception("FINAL GERMLINE ALIGNMENT END status=failed")
        raise
    finally:
        if query_path is not None:
            query_path.unlink(missing_ok=True)
        handler.close()
        audit.removeHandler(handler)


__all__ = [
    "FinalGermlineAlignmentError",
    "extract_final_sequences",
    "run_final_germline_alignment",
]
