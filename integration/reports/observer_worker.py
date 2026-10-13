"""Run the ORIGINAL back.py under a read-only Python execution observer.

No scientific rule is reimplemented. The observer records list changes after
the original assignment statements actually execute. Stdout stays untouched.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import runpy
import sys
from pathlib import Path

from integration.paths import BACK as CORE, require_virtual_environment
SUPPORTED_CORE_SHA256 = "da9a65c2bae69f162e9925ba13dcc5f672c94f354c80857248874accc16523d7"
RULE_IDS = ("murine_germline_difference", "cysteine_proline", "position_71", "fr4_motif")


def assignment_sites(path: Path) -> tuple[dict, int]:
    """Identify the six original write sites, not the conditions deciding mutations."""
    if hashlib.sha256(path.read_bytes()).hexdigest() != SUPPORTED_CORE_SHA256:
        raise ValueError("Observer does not support this back.py version; original file was not modified.")
    tree = ast.parse(path.read_text(encoding="utf-8"))
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "backmutate")
    parents = {child: parent for parent in ast.walk(function) for child in ast.iter_child_nodes(parent)}
    sites, initialization = {}, None
    for node in ast.walk(function):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if isinstance(target, ast.Name) and target.id == "backmutated":
            initialization = node.lineno
        if not isinstance(target, ast.Subscript) or not isinstance(target.value, ast.Name) or target.value.id != "backmutated":
            continue
        ancestors, current = [], node
        while current in parents:
            current = parents[current]
            ancestors.append(current)
        if any(isinstance(parent, ast.For) and isinstance(parent.iter, ast.Name) and parent.iter.id == "diffs" for parent in ancestors):
            rule = "murine_germline_difference"
        elif any(isinstance(parent, ast.If) and {"C", "P"}.issubset(
            {item.value for item in ast.walk(parent.test) if isinstance(item, ast.Constant) and isinstance(item.value, str)}
        ) for parent in ancestors):
            rule = "cysteine_proline"
        elif any(isinstance(parent, ast.If) and ast.unparse(parent.test) == "i == 71" for parent in ancestors):
            rule = "position_71"
        elif any(isinstance(part, ast.Name) and part.id == "FR4start" for part in ast.walk(target.slice)):
            rule = "fr4_motif"
        else:
            raise ValueError("Unrecognized original assignment; no rules will be fabricated.")
        sites[node.lineno] = rule
    expected = {"murine_germline_difference": 1, "cysteine_proline": 1, "position_71": 1, "fr4_motif": 3}
    if initialization is None or {rule: list(sites.values()).count(rule) for rule in RULE_IDS} != expected:
        raise ValueError("Original assignment sites do not match the supported core.")
    return sites, initialization


class ExecutionObserver:
    def __init__(self, path: Path):
        self.filename = str(path.resolve())
        self.sites, self.initialization_line = assignment_sites(path)
        self.candidates = []
        self.counters = {"VH": 0, "VL": 0}
        self.current = None
        self.record = None
        self.initialized = False
        self.await_initialization = False
        self.pending = None

    def _flush_assignment(self, frame):
        if self.pending is None:
            return
        candidate, rule, line, before = self.pending
        after = list(frame.f_locals["backmutated"])
        if len(before) != len(after):
            raise ValueError("Observer encountered a changing list length.")
        changes = [{"position": index, "from": old, "to": new}
                   for index, (old, new) in enumerate(zip(before, after), 1) if old != new]
        if len(changes) > 1:
            raise ValueError("One original scalar assignment changed multiple residues.")
        candidate["executed_assignments"][rule] += 1
        for change in changes:
            candidate["events"].append({"rule_id": rule, "core_line": line, **change})
        self.pending = None

    def trace(self, frame, event, arg):
        if frame.f_code.co_filename != self.filename or frame.f_code.co_name != "backmutate":
            return None
        if event == "call":
            self.current, self.record, self.pending = None, None, None
            self.initialized = self.await_initialization = False
            return self.trace
        if event not in {"line", "return"}:
            return self.trace
        self._flush_assignment(frame)
        local = frame.f_locals
        record = local.get("record")
        if record is not None and record is not self.record:
            chain = local["chain"]
            self.counters[chain] += 1
            self.record = record
            self.current = {"chain": chain, "candidate": self.counters[chain],
                            "proposal_id": record.id, "before_sequence": str(record.seq),
                            "after_sequence": None, "events": [],
                            "executed_assignments": {rule: 0 for rule in RULE_IDS}}
            self.candidates.append(self.current)
            self.initialized = self.await_initialization = False
        if self.await_initialization:
            self.initialized, self.await_initialization = True, False
        if self.current is not None and self.initialized:
            self.current["after_sequence"] = "".join(local["backmutated"])
        if event == "line" and self.current is not None:
            if frame.f_lineno == self.initialization_line:
                self.await_initialization = True
            if frame.f_lineno in self.sites and self.initialized:
                self.pending = (self.current, self.sites[frame.f_lineno], frame.f_lineno,
                                list(local["backmutated"]))
        return self.trace


def main() -> int:
    require_virtual_environment()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("vh_xml", type=Path)
    parser.add_argument("vl_xml", type=Path)
    parser.add_argument("observation_file", type=Path)
    args = parser.parse_args()
    observer = ExecutionObserver(CORE)
    original_argv, original_trace = sys.argv, sys.gettrace()
    try:
        sys.argv = [str(CORE), str(args.input_file), str(args.output_dir), str(args.vh_xml), str(args.vl_xml)]
        sys.settrace(observer.trace)
        runpy.run_path(str(CORE), run_name="__main__")
    finally:
        sys.settrace(original_trace)
        sys.argv = original_argv
    payload = {"core_file": str(CORE), "core_sha256": hashlib.sha256(CORE.read_bytes()).hexdigest(),
               "sites": {str(line): rule for line, rule in observer.sites.items()},
               "candidates": observer.candidates}
    args.observation_file.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
