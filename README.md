# python-bioinformatics

My first bioinformatics Python projects — a small collection of beginner scripts
for working with DNA, RNA, and protein sequences.

## Scripts

| Script | What it does |
| --- | --- |
| `base_counter.py` | Counts how many times each base (A, T, C, G) appears in a DNA sequence. |
| `DNA length calculator` | Returns the length (number of bases) of a DNA sequence. |
| `gc_halt.py` | Calculates the GC content of a DNA sequence, as a percentage. |
| `RNA Converter` | Transcribes DNA to RNA by replacing every T with U. |
| `Reverse DNA` | Reverses a DNA sequence. |
| `Reverse_complement(dna)` | Returns the reverse complement of a DNA sequence. |
| `rna_translation.py` | Translates an RNA sequence into a protein, showing 1-letter, 3-letter, and full amino acid names (stops at the first stop codon). |
| `orf_finder.py` | Notebook-friendly ORF finder that reads a FASTA file and finds open reading frames on both strands. Includes a minimum-length filter, maps reverse-strand coordinates back to the original sequence, collapses nested ORFs (one ORF per stop codon), and exports the results to a CSV file. |
| `Blast_annotate.py` | Add-on that takes the ORFs from `orf_finder.py` and BLASTs each predicted protein against NCBI to tell which are real, known genes and which are random ORFs. Requires `biopython` and an internet connection. |

## Running a script

Each script is self-contained and runs on its own. From the repository folder:

```bash
python base_counter.py
```

The example DNA/RNA sequence is defined at the bottom of each script — edit that
value to run it on your own sequence.

> **Note:** Some scripts print their output in Swedish (for example `DNA-sekvens`
> for "DNA sequence").

## ORF → gene annotation pipeline

`orf_finder.py` and `Blast_annotate.py` work together as a small genome
annotation pipeline: **FASTA → ORF prediction → BLAST annotation → CSV.**

1. **Find ORFs** (structural annotation) — locate the candidate genes.
2. **BLAST them** (functional annotation) — ask NCBI which candidates match a
   known protein, separating real genes from random ORFs.

First install Biopython (only needed for the BLAST step):

```bash
pip install biopython
```

Then, from a notebook or script in the repository folder:

```python
from orf_finder import run_orf_finder
from blast_annotate import annotate_orfs_with_blast, export_annotated_orfs_to_csv

# 1. Find ORFs in a genome
orfs = run_orf_finder("my_genome.fasta", min_protein_length=50)

# 2. BLAST the longest ORFs against NCBI (start small — each takes 1–3 minutes)
annotate_orfs_with_blast(orfs, max_orfs=3)

# 3. Save the annotated results
export_annotated_orfs_to_csv(orfs, "my_genome_annotation.csv")
```

The BLAST step prints a summary at the end, e.g.
`3 of 3 BLASTed ORFs matched a known protein (probable genes).`

> **Note:** online BLAST is slow and rate limited. Use `max_orfs` to test on a
> few ORFs first, and raise `min_protein_length` (e.g. to 100) on real genomes.
