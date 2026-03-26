# ViralTrace-MSA

ViralTrace-MSA is a modular toolkit for aligning closely related **circular viral genomes** using a **Center-Star progressive MSA** approach. It also computes per-column variability via a **PSSM** and **Shannon conservation entropy**.

## Directory structure

- `core/`
  - `genome_utils.py`: FASTA parsing and circular-genome calibration (rotation).
  - `alignment.py`: Global pairwise alignment and Center-Star MSA construction (including required gap propagation).
- `analysis/`
  - `statistics.py`: Build PSSM from an aligned FASTA and compute conservation entropy per alignment column.
- `tests/`
  - `test_genome_utils.py`: Unit tests for FASTA parsing, rotation, and calibration.
- `examples/`
  - `BLAST.fasta`: full genomes of some creatures
  - `example1.fasta`: small syntetic example
- `main.py`
  - CLI entry point running the pipeline: **Read -> Calibrate -> Align -> PSSM -> Plot**.

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Quick start (sample run)

Run on the bundled synthetic example:

```bash
python main.py \
  --input examples/BLAST.fasta \
  --output-alignment aligned_sample.fasta \
  --output-plot entropy_sample.png \
  --match 1 --mismatch -1 --gap -1 \
  --max-length 2000
```

Outputs:
- `aligned_sample.fasta`: the Center-Star MSA in FASTA format
- `entropy_sample.png`: conservation entropy (Shannon entropy) vs alignment position

## Running on your data (CLI)

```bash
python main.py \
  --input /path/to/input.fasta \
  --output-alignment /path/to/aligned.fasta \
  --output-plot /path/to/entropy.png \
  --match <match> --mismatch <mismatch> --gap <gap> \
  --max-length <N> \
  --random-seed <seed>
```

### Notes on runtime

Pairwise global alignment and calibration are computationally expensive for long sequences. Use `--max-length` to keep runs tractable.

## How the key steps work

### Circular genome calibration

`calibrate_genomes()` estimates a good rotation (start index) for each sequence using **local alignment** against a consensus. The goal is to orient genomes before MSA.

### Center-Star MSA

`perform_msa()`:
1. Computes all-pairs global alignment scores.
2. Chooses the center sequence with maximum total similarity score (C*).
3. Progressively aligns each sequence to the current center.
4. Propagates newly introduced **center gaps** into already aligned sequences so all aligned sequences end up with the same length.

### PSSM + conservation entropy

`build_pssm()` builds a frequency matrix (A/C/T/G/gap) per alignment column.
`compute_conservation_entropy()` computes Shannon entropy per column:
- low entropy = conserved columns
- high entropy = variable columns

## Testing

```bash
python -m pytest -q
```

