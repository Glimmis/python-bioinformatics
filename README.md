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
| `RNA -> Protein, RNA_translation.py` | Translates an RNA sequence into a protein, showing 1-letter, 3-letter, and full amino acid names (stops at the first stop codon). |
| `Open Reading Frame (ORF) Finder` | Finds open reading frames on both the forward strand and the reverse complement, and translates each ORF to protein. |
| `ORF Finder v2.0` | Extended ORF finder that also reports ORF length, GC content, protein length, and molecular weight, and exports the results to a CSV file. |

## Running a script

Each script is self-contained and runs on its own. From the repository folder:

```bash
python base_counter.py
```

The example DNA/RNA sequence is defined at the bottom of each script — edit that
value to run it on your own sequence.

> **Note:** Some scripts print their output in Swedish (for example `DNA-sekvens`
> for "DNA sequence").
