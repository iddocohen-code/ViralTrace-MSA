"""Statistical analysis utilities for aligned viral genomes."""

from __future__ import annotations

from math import log2

import numpy as np
import pandas as pd

from core.genome_utils import read_fasta_file

SequenceRecord = tuple[str, str]
ALPHABET = ["A", "C", "T", "G", "-"]


def build_pssm(input_path: str, alphabet: list[str] | None = None) -> pd.DataFrame:
    """Build a PSSM frequency matrix from an aligned FASTA file.

    Args:
        input_path: Path to an aligned FASTA file.
        alphabet: Optional ordered list of symbols to track.

    Returns:
        DataFrame with rows=symbols and columns=alignment positions.
    """
    symbols = alphabet if alphabet is not None else ALPHABET
    records = read_fasta_file(input_path)
    if not records:
        raise ValueError("PSSM cannot be built from an empty FASTA file.")

    aligned_sequences = [sequence for _, sequence in records]
    seq_lengths = {len(sequence) for sequence in aligned_sequences}
    if len(seq_lengths) != 1:
        raise ValueError("All aligned sequences must have identical lengths.")

    seq_len = len(aligned_sequences[0])
    num_sequences = len(aligned_sequences)
    matrix = np.zeros((len(symbols), seq_len), dtype=float)
    symbol_to_index = {symbol: index for index, symbol in enumerate(symbols)}

    for col in range(seq_len):
        counts = {symbol: 0 for symbol in symbols}
        for sequence in aligned_sequences:
            symbol = sequence[col]
            if symbol in counts:
                counts[symbol] += 1
        for symbol, count in counts.items():
            matrix[symbol_to_index[symbol], col] = count / num_sequences

    return pd.DataFrame(
        matrix,
        index=symbols,
        columns=[f"column_{i + 1}" for i in range(seq_len)],
    )


def compute_conservation_entropy(
    msa_data: list[SequenceRecord] | pd.DataFrame,
) -> pd.Series:
    """Compute per-column Shannon entropy for an MSA.

    Args:
        msa_data: Either MSA list of tuples `(id, aligned_sequence)` or a PSSM DataFrame.

    Returns:
        Series indexed by alignment columns with entropy values.
    """
    if isinstance(msa_data, pd.DataFrame):
        pssm = msa_data
    else:
        if not msa_data:
            return pd.Series(dtype=float)
        sequences = [sequence for _, sequence in msa_data]
        lengths = {len(sequence) for sequence in sequences}
        if len(lengths) != 1:
            raise ValueError("MSA sequences must have equal lengths.")
        pssm = _pssm_from_msa_records(msa_data)

    entropy_by_column: dict[str, float] = {}
    for column in pssm.columns:
        probabilities = pssm[column].to_numpy(dtype=float)
        entropy = 0.0
        for p in probabilities:
            if p > 0:
                entropy -= p * log2(p)
        entropy_by_column[column] = entropy
    return pd.Series(entropy_by_column, name="entropy")


def _pssm_from_msa_records(msa_records: list[SequenceRecord]) -> pd.DataFrame:
    """Build PSSM from in-memory MSA records."""
    sequences = [sequence for _, sequence in msa_records]
    seq_len = len(sequences[0])
    matrix = np.zeros((len(ALPHABET), seq_len), dtype=float)
    symbol_to_index = {symbol: index for index, symbol in enumerate(ALPHABET)}

    for col in range(seq_len):
        counts = {symbol: 0 for symbol in ALPHABET}
        for sequence in sequences:
            symbol = sequence[col]
            if symbol in counts:
                counts[symbol] += 1
        for symbol, count in counts.items():
            matrix[symbol_to_index[symbol], col] = count / len(sequences)

    return pd.DataFrame(
        matrix,
        index=ALPHABET,
        columns=[f"column_{i + 1}" for i in range(seq_len)],
    )

