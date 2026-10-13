"""Input parsing and sequential coordination shared by the provided dashboards."""
from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class SequencePair:
    name: str
    vh: str
    vl: str


def parse_sequence_pairs(content: str, validate: Callable) -> list[SequencePair]:
    content = content.lstrip("\ufeff").strip()
    if not content:
        raise ValueError("Batch input is empty.")
    first_line = content.splitlines()[0]
    delimiter = "\t" if "\t" in first_line else ";" if ";" in first_line else ","
    rows = list(csv.reader(io.StringIO(content), delimiter=delimiter, strict=True))
    pairs, names = [], set()
    for number, row in enumerate(rows, 1):
        if not row or not any(field.strip() for field in row):
            continue
        row = [field.strip() for field in row]
        if number == 1 and [field.lower() for field in row] == ["name", "vh", "vl"]:
            continue
        if len(row) != 3:
            raise ValueError(f"Row {number}: expected name, VH, VL.")
        name, vh, vl = row
        if not name or name in names:
            raise ValueError(f"Row {number}: sample name is empty or duplicated.")
        try:
            vh, vl = validate(vh, vl)
        except ValueError as exc:
            raise ValueError(f"Row {number} ({name}): {exc}") from exc
        names.add(name)
        pairs.append(SequencePair(name, vh, vl))
    if not pairs:
        raise ValueError("No sequence pairs found in batch input.")
    return pairs


def run_sequentially(run_ids, run_one: Callable) -> None:
    for run_id in run_ids:
        run_one(run_id)
