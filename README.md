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
  --input examples/example1.fasta \
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

### Circular Genome Calibration (Heuristic Algorithm)

Since viral genomes are often circular, their starting positions in FASTA files are arbitrary. To ensure a meaningful MSA, the toolkit implements a sophisticated iterative heuristic to co-orient the sequences:

1.  **Initial Consensus:** A primary consensus sequence is generated from the input records to serve as the initial reference.
2.  **Iterative Re-orientation:**
    * **Local Alignment:** Each genome is locally aligned against the current reference to identify the region of maximal similarity.
    * **Rotation:** Each sequence is rotated so that the start of the highest-scoring local alignment becomes the new sequence start.
    * **Evaluation:** The algorithm calculates the total similarity score of the newly rotated set.
3.  **Adaptive Optimization:**
    * **Greedy Step:** If the total score improves significantly, the new orientation is accepted, and a new consensus is generated.
    * **Stochastic Escape:** If no improvement is found, a "Certainty Tracker" (stagnation counter) increases. A random sequence from the set is then chosen as the new reference to "shake" the algorithm out of local optima and discover better global alignments.
4.  **Termination:** The process concludes when the orientation remains stable for a threshold of iterations (10x the number of sequences), ensuring a high-confidence starting point for the subsequent MSA.
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

