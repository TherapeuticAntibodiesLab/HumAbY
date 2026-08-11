import pytest

from batch_analysis import parse_sequence_pairs, run_sequentially


VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARWGQGTLVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYWSTPITFGQGTKVEIK"


def validator(vh, vl):
    if len(vh) < 50 or len(vl) < 50:
        raise ValueError("too short")
    return vh.upper(), vl.upper()


def test_parse_csv_and_tsv_with_optional_header():
    csv_pairs = parse_sequence_pairs(f"name,VH,VL\na,{VH},{VL}\nb,{VH},{VL}\n", validator)
    tsv_pairs = parse_sequence_pairs(f"a\t{VH}\t{VL}\n", validator)
    assert [pair.name for pair in csv_pairs] == ["a", "b"]
    assert tsv_pairs[0].vh == VH


def test_duplicate_names_are_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        parse_sequence_pairs(f"a,{VH},{VL}\nA,{VH},{VL}\n", validator)


def test_runner_is_strictly_sequential():
    events = []
    run_sequentially(["one", "two"], lambda item: events.extend([f"start:{item}", f"end:{item}"]))
    assert events == ["start:one", "end:one", "start:two", "end:two"]
