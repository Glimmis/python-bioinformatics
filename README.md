# python-bioinformatics

A small genome annotation pipeline in Python: **FASTA → ORF prediction → BLAST
annotation → CSV.**

Given a genome or contig, it finds the candidate protein-coding genes and then
asks NCBI which of them correspond to proteins that are actually known — the
two halves of gene annotation, structural and functional, in about 700 lines of
dependency-light Python.

```bash
# 1. Structural annotation: find the ORFs
python orf_finder.py examples/example_contig.fasta --min-length 100

# 2. Functional annotation: BLAST them against NCBI nr
python blast_annotate.py example_contig_orfs.csv --max-orfs 3
```

## Why

An ORF finder on its own reports every stretch between a start and a stop
codon, and most of those occur by chance. The interesting question is which
candidates are real. Running the predicted proteins through BLASTp against the
non-redundant database separates the two, so the output is a ranked list of
probable genes rather than a list of possibilities.

## The two steps

### `orf_finder.py` — structural annotation

Scans every sequence in a FASTA file for complete ORFs (ATG to stop codon) in
all three reading frames on both strands, translates them, and writes a CSV
with coordinates, strand, frame, length, GC content, protein sequence and
estimated molecular weight.

Details worth knowing about:

- **Reverse-strand coordinates are mapped back to the original sequence**, so
  start/end always refer to the numbering in the input FASTA rather than to the
  reverse-complemented copy.
- **ORFs sharing a stop codon are collapsed to the longest one.** Such ORFs are
  the same gene read from different internal start codons; reporting all of
  them would inflate the count several-fold. Real ORF finders report one ORF
  per stop codon, and so does this one.
- **The length filter is applied before translation**, so large inputs do not
  pay to translate ORFs that will then be discarded.
- Reverse complement uses `str.translate` — O(n) rather than the O(n²) of
  character-by-character concatenation.

### `blast_annotate.py` — functional annotation

Takes the ORF CSV and runs BLASTp against NCBI `nr` for each predicted protein,
longest first. It records the hit title, E-value and percent identity, and
classifies each ORF as a probable gene or probable noise.

Hits titled "hypothetical", "uncharacterized" or "unnamed protein" are skipped
in favour of a properly named hit further down the list, with the best
significant hit retained as a fallback so a real gene is never lost to a
placeholder annotation.

Online BLAST is slow (roughly 0.5–3 minutes per sequence) and rate limited, so
`--max-orfs` is there to test on a handful first. For hundreds of sequences,
install BLAST+ locally instead.

## Worked example

`examples/example_contig.fasta` is a synthetic 3 612 bp contig with four
protein-coding genes planted at known positions — three on the forward strand,
one on the reverse — inside random intergenic sequence. It is generated from a
fixed random seed by `examples/make_example_contig.py`, so it is reproducible.

```bash
python orf_finder.py examples/example_contig.fasta --min-length 100 \
    -o examples/example_contig_orfs.csv
```

The committed output is `examples/example_contig_orfs.csv`:

| strand | frame | start | end | bp | aa | GC % | longest |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Forward | 1 | 994 | 1716 | 723 | 240 | 48.27 | |
| Forward | 2 | 251 | 793 | 543 | 180 | 46.78 | |
| Forward | 2 | 2480 | 3412 | 933 | 310 | 53.38 | ✓ |
| Reverse complement | 2 | 1917 | 2279 | 363 | 120 | 51.24 | |

All four planted genes are recovered at the right coordinates and on the right
strand, and the nine raw ORFs collapse to four. The reverse-strand gene is
reported in forward-sequence coordinates, which is the coordinate mapping doing
its job.

The proteins are synthetic, so BLASTing them returns nothing — the annotation
step is meant to be run on real sequence.

## Real data: phiX174

`examples/phix174.fasta` is the real phiX174 genome (NC_001422.1, 5,386 bp),
downloaded from NCBI. Unlike the synthetic contig above, its proteins are real
and BLASTing them returns actual matches.

```bash
python orf_finder.py examples/phix174.fasta --min-length 50 \
    -o examples/phix174_orfs.csv
python blast_annotate.py examples/phix174_orfs.csv --max-orfs 5 \
    -o examples/phix174_annotation.csv
```

`orf_finder.py` found 22 ORFs (16 forward, 6 reverse complement) at
`--min-length 50`; the five longest were BLASTed. The committed output is
`examples/phix174_annotation.csv`:

| strand | start | end | aa | BLAST hit | E-value |
| --- | --- | --- | --- | --- | --- |
| Forward | 1001 | 2284 | 427 | MULTISPECIES: major capsid protein [Bacteria] (matches capsid protein F, Escherichia phage phiX174) | 0.0 |
| Reverse complement | 1094 | 1651 | 185 | conserved hypothetical protein [Mesorhizobium sp. ORS 3324] | 1.92504e-128 |
| Forward | 2395 | 2922 | 175 | No significant hit | |
| Forward | 2931 | 3917 | 328 | MULTISPECIES: Minor spike protein [Bacteria] (matches minor spike protein H, Sinsheimervirus phiX174) | 0.0 |
| Forward | 3076 | 3684 | 202 | Minor spike protein H (modular protein) [Microbacterium sp. C448] | 2.7869e-110 |

Two of the five BLASTed ORFs — the 427 aa and 328 aa forward-strand proteins —
matched known phiX174 genes by name (capsid protein F and minor spike protein
H); the other two significant hits were homologs in unrelated organisms, and
the 175 aa ORF returned no significant hit.

## Installation

The ORF finder needs only the standard library. The BLAST step needs Biopython:

```bash
pip install -r requirements.txt
```

## Command-line options

```
orf_finder.py FASTA [-o CSV] [-m MIN_LENGTH] [--detail] [--no-csv] [-q]
blast_annotate.py ORF_CSV [-o CSV] [-m MIN_LENGTH] [-n MAX_ORFS] [-e EVALUE] [--pause SECONDS]
```

Both modules are importable as well:

```python
from orf_finder import run_orf_finder
from blast_annotate import annotate_orfs_with_blast, export_annotated_orfs_to_csv

orfs = run_orf_finder("genome.fasta", min_protein_length=100)
annotate_orfs_with_blast(orfs, max_orfs=3)
export_annotated_orfs_to_csv(orfs, "genome_annotation.csv")
```

## Repository layout

```
orf_finder.py           structural annotation
blast_annotate.py       functional annotation
examples/               synthetic contig, its ORF output, and the generator script
basics/                 short single-purpose scripts written while learning Python
```

`basics/` holds the early exercises — base counting, GC content, reverse
complement, translation — kept for reference.
