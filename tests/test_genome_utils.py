"""Unit tests for genome utility logic."""

import random

from core.genome_utils import calibrate_genomes, read_fasta_file, rotate_seq_from_point


def test_read_fasta_file_parses_records(tmp_path):
    fasta_path = tmp_path / "toy.fasta"
    fasta_path.write_text(">s1\nACTG\n>s2\nA-TA\n", encoding="utf-8")

    result = read_fasta_file(str(fasta_path))
    assert result == [("s1", "ACTG"), ("s2", "A-TA")]


def test_rotate_seq_from_point_wraps():
    assert rotate_seq_from_point("ACTG", 2) == "TGAC"
    assert rotate_seq_from_point("ACTG", 6) == "TGAC"


def test_calibrate_genomes_noop_for_single_sequence():
    records = [("s1", "ACTGACTG")]
    assert calibrate_genomes(records) == records


def test_calibrate_genomes_returns_same_ids_and_lengths():
    records = [
        ("a", "ACTGACTG"),
        ("b", "TGACTGAC"),
        ("c", "GACTGACT"),
    ]
    calibrated = calibrate_genomes(records, rng=random.Random(123))

    assert [record_id for record_id, _ in calibrated] == ["a", "b", "c"]
    assert [len(seq) for _, seq in calibrated] == [8, 8, 8]

