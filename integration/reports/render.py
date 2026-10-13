"""The single official compact TXT representation of complete report data."""
from __future__ import annotations

from .reader import RULES_UNAVAILABLE

_OMITTED_LINES = frozenset({
    "Organized from existing outputs; sequences and original reports are unchanged.",
    "Murine database version/source recorded in run metadata: Not recorded for this execution",
    "Reference labels below distinguish the human IMGT grafting template from the murine NCBI/IMGT Backmutation reference.",
    "Historical origin is only named when verified; record codes are not automatically GenBank accessions.",
    "Final backmt germline annotations below are ADDITIONAL BLAST matches to the human database, not original template IDs or definitive assignments.",
    "Database origin: IMGT-derived human reconstructed frameworks",
    "Rule counts are recorded during the verified replay; no-op writes have zero changes and rule positions may overlap.",
})
_OMITTED_PREFIXES = (
    "Identifier type:", "Reference evidence:", "Interpretation:",
    "Ties:", "Analysis scope:", "Rule evidence:",
)


def compact_existing_text(text: str) -> str:
    """Use this policy only to normalize already-saved TXT, never as scientific input."""
    return "".join(line for line in text.splitlines(keepends=True)
                   if line.strip() not in _OMITTED_LINES
                   and not line.strip().startswith(_OMITTED_PREFIXES))


def render_text(document: dict) -> str:
    lines = ["BACKMUTATION RESULTS", "====================", "",
             "Positions: " + document["position_numbering"],
             "Run: " + str(document["run"]["run_id"] or "Not recorded"),
             "Human database path recorded for this execution: " + str(document["run"]["database_path"] or "Not recorded")]
    murine = document["run"].get("murine_database_source")
    if murine and murine != "Not recorded for this execution":
        lines.append("Murine database version/source recorded in run metadata: " + str(murine))
    marks = document.get("comparison")
    if marks:
        lines += ["Result SHA-256: " + marks["result_sha256"],
                  "Execution SHA-256: " + (marks["execution_sha256"] or "Not recorded for this execution")]
    lines += [f"Total candidates: {document['totals']['candidates']}",
              f"Total effective back-mutations: {document['totals']['effective_changes']}", ""]
    for name, chain in document["chains"].items():
        reference = chain.get("mouse_reference", {})
        origin = reference.get("origin") or "origin not verified"
        lines += ["CHAIN " + name, "========", "", "Framework map:",
                  chain["framework_map"] or "  Not available", "",
                  f"Murine reference for Backmutation (mGerm) — {origin}:",
                  "  Sequence reference code: " + (reference.get("identifier") or "Not available"),
                  "  Reference accession field: " + (reference.get("accession") or "Not available"),
                  "  Database path saved in XML: " + (reference.get("database_path_recorded_in_xml") or "Not available"),
                  "  mGerm: " + (chain["mouse_germline"] or "Not available"),
                  "  match: " + (chain["match"] or "Not available"), "", "Input sequence:",
                  "  " + (chain["input_label"] or "input_") + ": " + (chain["input_sequence"] or "Not available"), ""]
        for candidate in chain["candidates"]:
            title = f"{name} candidate {candidate['number']}"
            reference = candidate.get("human_reference", {})
            origin = reference.get("origin") or "origin not verified"
            genes = " ".join(f"{key}:{value}" for key, value in reference.get("genes", {}).items())
            lines += [title, "-" * len(title), f"Human reference used for grafting (hGerm) — {origin}:",
                      "  Framework reference code: " + (reference.get("identifier") or "Not available"),
                      "  Saved candidate FASTA ID: " + reference.get("saved_fasta_id", candidate["proposal_id"]),
                      "  IMGT segment names: " + (genes or "Not available"),
                      "  hGerm: " + (candidate["human_germline"] or "Not available"),
                      "Humanized sequence before Backmutation:", "  proposed: " + candidate["humanized_sequence"],
                      "Backmutated sequence:", "  backmt: " + candidate["backmutated_sequence"],
                      f"Back-mutations performed: {candidate['count']}",
                      "Changes (position: humanized -> backmutated):"]
            lines += [f"  {change['position']}: {change['from']} -> {change['to']}" for change in candidate["changes"]]
            if not candidate["changes"]:
                lines.append("  No effective residue changes.")
            lines += ["", "Germline determined from this backmt sequence:"]
            annotation = candidate.get("germline_annotation")
            if annotation:
                lines += [f"  {key}: {value}" for key, value in annotation.items()
                          if key not in {"Database origin", "Interpretation", "Ties", "Analysis scope"}]
            else:
                lines.append("  Unavailable; no additional analysis was performed.")
            lines += ["", "Rules applied:"]
            if candidate["rules_available"]:
                lines.append("  Rule | Changes | Positions")
                lines += [f"  {rule['name']} | {rule['count']} | {rule['positions']}" for rule in candidate["rules"]]
            else:
                lines.append("  Rule evidence unavailable: " + candidate.get("rules_error", "original logs and frozen observation inputs do not provide validated rules."))
                lines.append("  " + RULES_UNAVAILABLE)
            lines.append("")
    # Source hashes and detailed analysis provenance are kept in JSON, not repeated in TXT.
    analysis = document.get("final_germline_analysis", {})
    if analysis.get("error"):
        lines += ["FINAL GERMLINE ANALYSIS", "======================", "Error: " + analysis["error"]]
    return "\n".join(lines) + "\n"
