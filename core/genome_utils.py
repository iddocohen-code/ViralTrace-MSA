"""Genome utility functions for ViralTrace-MSA."""

from __future__ import annotations

import random
from collections import Counter
from itertools import zip_longest
from typing import Optional

from Bio.Align import PairwiseAligner

SequenceRecord = tuple[str, str]


def read_fasta_file(path: str) -> list[SequenceRecord]:
    """Read sequences from a FASTA file.

    Args:
        path: Path to a FASTA file.

    Returns:
        A list of tuples in the form `(sequence_id, sequence_string)`.
    """
    records: list[SequenceRecord] = []
    with open(path, "r", encoding="utf-8") as handle:
        current_id: Optional[str] = None
        current_sequence_parts: list[str] = []

        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_sequence_parts)))
                current_id = line[1:].strip()
                current_sequence_parts = []
            else:
                current_sequence_parts.append(line)

        if current_id is not None:
            records.append((current_id, "".join(current_sequence_parts)))
    return records


def new_consensus_seq(
    sequences: list[str], rng: Optional[random.Random] = None
) -> str:
    """Build a consensus sequence using column voting.

    Args:
        sequences: Input sequences for consensus creation.
        rng: Optional random generator used for tie-breaking.

    Returns:
        Consensus sequence string.
    """
    if not sequences:
        return ""

    randomizer = rng if rng is not None else random
    consensus: list[str] = []

    for column in zip_longest(*sequences, fillvalue="N"):
        counts = Counter(column)
        counts.pop("N", None)
        if not counts:
            consensus.append("N")
        elif len(set(counts.values())) == 1 and len(counts) > 1:
            consensus.append(randomizer.choice(list(counts.keys())))
        else:
            consensus.append(counts.most_common(1)[0][0])

    return "".join(consensus)


def find_starting_point_by_local_alignment(
    reference_seq: str, current_seq: str
) -> tuple[int, float]:
    """Find the best local-alignment start point for circular calibration.

    Args:
        reference_seq: Reference sequence.
        current_seq: Sequence to align against the reference.

    Returns:
        Tuple `(start_point, score)` from the best local alignment.
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1

    best_alignment = aligner.align(reference_seq, current_seq)[0]
    similarity_score = float(best_alignment.score)
    start_point = int(best_alignment.aligned[1][0][0])
    return start_point, similarity_score


def rotate_seq_from_point(sequence: str, point: int) -> str:
    """Rotate a sequence so `point` becomes the first position."""
    if not sequence:
        return sequence
    normalized = point % len(sequence)
    return sequence[normalized:] + sequence[:normalized]


def calibrate_genomes(
    seq_record_list: list[SequenceRecord],
    min_improvement: float = 0.001,
    max_stagnation_factor: int = 10,
    rng: Optional[random.Random] = None,
) -> list[SequenceRecord]:
    """Calibrate circular genomes by iterative local-alignment-guided rotation.

    Args:
        seq_record_list: Input sequence records.
        min_improvement: Minimum average score gain to accept an iteration.
        max_stagnation_factor: Max stagnation multiplier of sequence count.
        rng: Optional random generator for deterministic fallback/reference choice.

    Returns:
        Calibrated sequence records in the same `(id, sequence)` format.
    """
    if len(seq_record_list) < 2:
        return seq_record_list

    randomizer = rng if rng is not None else random
    best_sequences = list(seq_record_list)
    reference_seq = new_consensus_seq([seq for _, seq in best_sequences], rng=randomizer)

    stagnation_counter = 0
    max_score = float("-inf")
    max_stagnation = len(best_sequences) * max_stagnation_factor

    while stagnation_counter < max_stagnation:
        rotated_sequences: list[SequenceRecord] = []
        total_score = 0.0

        for sequence_id, sequence in best_sequences:
            start_point, alignment_score = find_starting_point_by_local_alignment(
                reference_seq, sequence
            )
            rotated = rotate_seq_from_point(sequence, start_point)
            rotated_sequences.append((sequence_id, rotated))
            total_score += alignment_score

        average_score = total_score / len(best_sequences)
        if average_score > max_score + min_improvement:
            max_score = average_score
            stagnation_counter = 0
            best_sequences = rotated_sequences
            reference_seq = new_consensus_seq(
                [seq for _, seq in best_sequences], rng=randomizer
            )
        else:
            stagnation_counter += 1
            reference_seq = randomizer.choice(best_sequences)[1]

    return best_sequences

