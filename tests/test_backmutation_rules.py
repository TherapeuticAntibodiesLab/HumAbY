from dash_app import parse_backmutation_display
from humanize import enrich_backmutation_report


def test_web_parser_reads_compact_rule_table():
    report = (
        "VH candidate 1\n"
        "Back-mutations performed: 2\n"
        "Rules applied:\n"
        "  Rule | Changes | Positions\n"
        "  Cysteine/Proline protection | 1 | 48\n"
        "  FR4 motif correction | 1 | 104\n"
    )
    parsed = parse_backmutation_display(report)[("VH", 1)]
    assert parsed["rules"] == [
        {"name": "Cysteine/Proline protection", "count": 1, "positions": "48"},
        {"name": "FR4 motif correction", "count": 1, "positions": "104"},
    ]


def test_text_report_organizes_machine_rules(monkeypatch):
    raw = (
        "Chain VH\n"
        "hGerm HUMAN\n"
        "backmt ACDEFGHIKLMNPQRSTVWY\n"
        "rule murine_germline_difference count=2 positions=12,48\n"
        "rule cysteine_proline count=0 positions=-\n"
        "rule position_71 count=1 positions=72\n"
        "rule fr4_motif count=0 positions=-\n"
    )
    monkeypatch.setattr("humanize._format_backmutated_germline_info", lambda *args: ["  germline: test"])
    report = enrich_backmutation_report(raw, "human-db", proposed_sequences={"VH": [], "VL": []})
    assert report.count("Rules applied:") == 1
    assert "Murine vs murine germline difference | 2 | 12,48" in report
    assert "Internal position 71 restoration | 1 | 72" in report
