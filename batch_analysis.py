"""Independent helpers for sequential paired-sequence batch analysis."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from typing import Callable, Iterable, List


@dataclass(frozen=True)
class SequencePair:
    name: str
    vh: str
    vl: str


def parse_sequence_pairs(
    text: str,
    validator: Callable[[str, str], tuple[str, str]],
) -> List[SequencePair]:
    """Parse CSV/TSV rows formatted as name,VH,VL and validate every pair."""
    if not text or not text.strip():
        raise ValueError("The batch input is empty.")
    sample = text[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
    except csv.Error:
        dialect = csv.excel
    rows = list(csv.reader(io.StringIO(text), dialect))
    rows = [row for row in rows if row and any(cell.strip() for cell in row)]
    if not rows:
        raise ValueError("The batch input contains no sequence pairs.")

    first = [cell.strip().lower() for cell in rows[0]]
    if len(first) >= 3 and first[0] in {"name", "nombre", "id"} and first[1:3] == ["vh", "vl"]:
        rows = rows[1:]
    if not rows:
        raise ValueError("The batch input has a header but no sequence pairs.")

    pairs = []
    seen_names = set()
    for number, row in enumerate(rows, 1):
        if len(row) != 3:
            raise ValueError(f"Batch row {number} must contain exactly: name, VH, VL.")
        name, vh_text, vl_text = (cell.strip() for cell in row)
        if not name:
            raise ValueError(f"Batch row {number} has no name.")
        normalized_name = name.casefold()
        if normalized_name in seen_names:
            raise ValueError(f"Duplicate batch name: {name}.")
        try:
            vh, vl = validator(vh_text, vl_text)
        except ValueError as exc:
            raise ValueError(f"Invalid pair '{name}': {exc}") from exc
        seen_names.add(normalized_name)
        pairs.append(SequencePair(name=name, vh=vh, vl=vl))
    return pairs


def run_sequentially(items: Iterable[str], runner: Callable[[str], None]) -> None:
    """Run each identifier synchronously and in input order."""
    for item in items:
        runner(item)


__all__ = ["SequencePair", "parse_sequence_pairs", "run_sequentially"]
