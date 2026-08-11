import subprocess

import pytest

from final_germline_alignment import (
    FinalGermlineAlignmentError,
    extract_final_sequences,
    run_final_germline_alignment,
)


RAW_OUTPUT = """
 Chain VH
backmt EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARWGQGTLVTVSS

 Chain VL
backmt DIQMTQSPSSLSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYWSTPITFGQGTKVEIK
"""


def test_extract_final_sequences_tracks_chain_and_candidate_number():
    sequences = extract_final_sequences(RAW_OUTPUT)
    assert [item[0] for item in sequences] == ["VH_final_1", "VL_final_1"]
    assert all(set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY") for _, sequence in sequences)


def test_extract_final_sequences_requires_chain_header():
    with pytest.raises(FinalGermlineAlignmentError):
        extract_final_sequences("backmt ACDEFGHIKLMNPQRSTVWY")


def test_run_writes_top_five_report_and_audit_log(tmp_path, monkeypatch):
    sequence = extract_final_sequences(RAW_OUTPUT)[0][1]
    rows = []
    for rank in range(1, 6):
        rows.append(
            "\t".join([
                "VH_final_1", f"VH_{rank:06d}", "99.000", str(len(sequence)),
                "1", "0", "1", str(len(sequence)), "1", str(len(sequence)),
                f"1e-{50-rank}", str(300-rank), sequence, sequence,
                f"Heavy chain germline V:IGHV{rank}-1*01 J:IGHJ4*02",
            ])
        )

    def fake_run(cmd, **kwargs):
        assert cmd[0] == "blastp"
        assert cmd[cmd.index("-max_target_seqs") + 1] == "5"
        return subprocess.CompletedProcess(cmd, 0, "\n".join(rows), "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    report = run_final_germline_alignment(RAW_OUTPUT, "human_db", str(tmp_path))

    text = report.read_text(encoding="utf-8")
    assert "QUERY: VH_final_1" in text
    assert "QUERY: VL_final_1" in text
    assert "#5  VH_000005" in text
    assert "Query     " in text and "Germline  " in text
    log = (tmp_path / "final_germline_alignment.log").read_text(encoding="utf-8")
    assert "status=success" in log
