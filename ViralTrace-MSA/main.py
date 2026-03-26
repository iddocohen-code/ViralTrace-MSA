"""CLI entry point for ViralTrace-MSA."""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from analysis.statistics import build_pssm, compute_conservation_entropy
from core.alignment import perform_msa
from core.genome_utils import calibrate_genomes, read_fasta_file

SequenceRecord = tuple[str, str]


def write_msa_to_fasta(msa: list[SequenceRecord], path: str) -> None:
    """Write an MSA list of tuples to a FASTA file.

    Args:
        msa: Aligned sequence records as `(sequence_id, aligned_sequence)`.
        path: Destination FASTA path.
    """
    records = [
        SeqRecord(Seq(sequence), id=sequence_id, description="")
        for sequence_id, sequence in msa
    ]
    alignment = MultipleSeqAlignment(records)
    AlignIO.write(alignment, path, "fasta")


def run_pipeline(
    input_path: str,
    output_alignment_path: str,
    output_plot_path: str | None,
    match: float,
    mismatch: float,
    gap: float,
    max_length: int,
    random_seed: int | None,
) -> tuple[list[SequenceRecord], pd.DataFrame, pd.Series]:
    """Run the full ViralTrace-MSA workflow.

    Pipeline: Read -> Calibrate -> Align -> PSSM -> Plot.

    Args:
        input_path: Input FASTA file path.
        output_alignment_path: Output FASTA path for aligned sequences.
        output_plot_path: Optional output path for entropy plot.
        match: Pairwise alignment match score.
        mismatch: Pairwise alignment mismatch penalty.
        gap: Pairwise alignment gap penalty (open == extend).
        max_length: Trim sequences to this length before calibration+MSA.
        random_seed: Optional seed for deterministic calibration tie-breaking.

    Returns:
        A tuple `(msa_records, pssm_df, entropy_series)`.
    """
    rng = None
    if random_seed is not None:
        import random

        rng = random.Random(random_seed)

    # Important performance choice:
    # Calibration uses expensive local alignment, so we bound runtime by trimming
    # before calibration.
    records = read_fasta_file(input_path)
    trimmed_records = [(sid, seq[:max_length]) for sid, seq in records]

    calibrated = calibrate_genomes(trimmed_records, rng=rng)
    msa = perform_msa(calibrated, match=match, mismatch=mismatch, gap=gap)
    write_msa_to_fasta(msa, output_alignment_path)

    pssm = build_pssm(output_alignment_path)
    entropy = compute_conservation_entropy(pssm)

    if output_plot_path:
        _plot_entropy(entropy, output_plot_path)

    return msa, pssm, entropy


def _plot_entropy(entropy: pd.Series, output_plot_path: str) -> None:
    """Create and save an entropy profile plot."""
    plot_df = pd.DataFrame(
        {"position": range(1, len(entropy) + 1), "entropy": entropy.to_numpy()}
    )
    plt.figure(figsize=(12, 4))
    sns.lineplot(data=plot_df, x="position", y="entropy")
    plt.title("Per-position Conservation Entropy")
    plt.xlabel("Alignment Position")
    plt.ylabel("Shannon Entropy")
    plt.tight_layout()
    plt.savefig(output_plot_path, dpi=200)
    plt.close()


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="ViralTrace-MSA: calibration, Center-Star alignment, and analysis."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file path.")
    parser.add_argument(
        "--output-alignment",
        required=True,
        help="Output FASTA path for aligned sequences.",
    )
    parser.add_argument(
        "--output-plot",
        default=None,
        help="Optional output path for entropy plot (e.g. entropy.png).",
    )
    parser.add_argument("--match", type=float, default=1.0, help="Match score.")
    parser.add_argument(
        "--mismatch", type=float, default=-1.0, help="Mismatch penalty."
    )
    parser.add_argument("--gap", type=float, default=-1.0, help="Gap penalty.")
    parser.add_argument(
        "--max-length",
        type=int,
        default=10000,
        help="Trim sequences to this length before calibration+MSA.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="Optional seed for deterministic calibration behavior.",
    )
    return parser


def main() -> None:
    """Run CLI."""
    parser = build_arg_parser()
    args = parser.parse_args()
    run_pipeline(
        input_path=args.input,
        output_alignment_path=args.output_alignment,
        output_plot_path=args.output_plot,
        match=args.match,
        mismatch=args.mismatch,
        gap=args.gap,
        max_length=args.max_length,
        random_seed=args.random_seed,
    )


if __name__ == "__main__":
    main()

