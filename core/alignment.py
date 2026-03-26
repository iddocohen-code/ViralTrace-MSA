"""Center-Star MSA alignment logic."""

from __future__ import annotations

from itertools import combinations

from Bio import pairwise2

SequenceRecord = tuple[str, str]
PairwiseAlignment = tuple[str, str, float]


def pairwise_align(
    seq1: SequenceRecord,
    seq2: SequenceRecord,
    match: float,
    mismatch: float,
    gap: float,
) -> PairwiseAlignment:
    """Perform global pairwise alignment and return the top-scoring result.

    Args:
        seq1: First sequence tuple `(id, sequence)`.
        seq2: Second sequence tuple `(id, sequence)`.
        match: Match score.
        mismatch: Mismatch score.
        gap: Gap open and extension penalty.

    Returns:
        A tuple `(aligned_seq1, aligned_seq2, score)`.
    """
    _, s1 = seq1
    _, s2 = seq2
    aligned_s1, aligned_s2, score, _, _ = pairwise2.align.globalms(
        s1, s2, match, mismatch, gap, gap, one_alignment_only=True
    )[0]
    return aligned_s1, aligned_s2, float(score)


def _all_pairwise_global_alignments(
    seq_record_list: list[SequenceRecord],
    match: float,
    mismatch: float,
    gap: float,
) -> dict[tuple[str, str], PairwiseAlignment]:
    """Compute global alignment for each unique pair."""
    results: dict[tuple[str, str], PairwiseAlignment] = {}
    for seq1, seq2 in combinations(seq_record_list, 2):
        results[(seq1[0], seq2[0])] = pairwise_align(
            seq1, seq2, match, mismatch, gap
        )
    return results


def _find_c_star(
    all_global_alignments: dict[tuple[str, str], PairwiseAlignment],
) -> SequenceRecord:
    """Find Center-Star sequence by maximum total pairwise score."""
    score_sums: dict[str, float] = {}
    ungapped_sequences: dict[str, str] = {}

    for (id1, id2), (aligned1, aligned2, score) in all_global_alignments.items():
        if id1 not in score_sums:
            score_sums[id1] = 0.0
            ungapped_sequences[id1] = aligned1.replace("-", "")
        if id2 not in score_sums:
            score_sums[id2] = 0.0
            ungapped_sequences[id2] = aligned2.replace("-", "")
        score_sums[id1] += score
        score_sums[id2] += score

    center_id = max(score_sums, key=score_sums.get)
    return center_id, ungapped_sequences[center_id]


def _merge_into_msa(
    msa: dict[str, str],
    center_id: str,
    aligned_center: str,
    aligned_sequence: str,
    sequence_id: str,
) -> None:
    """Merge a new center-vs-sequence alignment into existing MSA in-place.

    This performs the required Center-Star gap propagation:
    when the center alignment introduces gaps, those gaps are inserted into all
    previously aligned sequences.
    """
    old_center = msa[center_id]
    gap_positions: list[int] = []
    old_index = 0

    for new_index, char in enumerate(aligned_center):
        if old_index < len(old_center) and old_center[old_index] == char:
            old_index += 1
        elif char == "-":
            gap_positions.append(new_index)
        else:
            old_index += 1

    if gap_positions:
        for existing_id in list(msa.keys()):
            if existing_id == center_id:
                continue
            seq_list = list(msa[existing_id])
            for pos in sorted(gap_positions, reverse=True):
                seq_list.insert(pos, "-")
            msa[existing_id] = "".join(seq_list)

    msa[center_id] = aligned_center
    msa[sequence_id] = aligned_sequence


def _msa_builder(
    seq_record_list: list[SequenceRecord],
    c_star: SequenceRecord,
    match: float,
    mismatch: float,
    gap: float,
) -> list[SequenceRecord]:
    """Build progressive star MSA around chosen center sequence."""
    center_id, _ = c_star
    msa: dict[str, str] = {center_id: c_star[1]}

    for sequence_id, sequence in seq_record_list:
        if sequence_id == center_id:
            continue
        aligned_center, aligned_curr, _ = pairwise_align(
            (center_id, msa[center_id]),
            (sequence_id, sequence),
            match,
            mismatch,
            gap,
        )
        _merge_into_msa(
            msa=msa,
            center_id=center_id,
            aligned_center=aligned_center,
            aligned_sequence=aligned_curr,
            sequence_id=sequence_id,
        )

    return list(msa.items())


def perform_msa(
    seq_record_list: list[SequenceRecord],
    match: float,
    mismatch: float,
    gap: float,
) -> list[SequenceRecord]:
    """Perform Center-Star MSA on a list of sequence records.

    Args:
        seq_record_list: Sequence records `(id, sequence)` to align.
        match: Match score.
        mismatch: Mismatch score.
        gap: Gap open and extension penalty.

    Returns:
        MSA as list of `(id, aligned_sequence)` tuples.
    """
    if not seq_record_list:
        return []
    if len(seq_record_list) == 1:
        return seq_record_list

    alignments = _all_pairwise_global_alignments(seq_record_list, match, mismatch, gap)
    c_star = _find_c_star(alignments)
    return _msa_builder(seq_record_list, c_star, match, mismatch, gap)

