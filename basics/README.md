# basics

Short single-purpose scripts written while learning Python for sequence work.
They are kept here for reference; the working code in this repository is the
ORF-to-annotation pipeline in the parent directory.

| Script | What it does |
| --- | --- |
| `base_counter.py` | Counts A, T, C and G in a DNA sequence. |
| `dna_length.py` | Returns the length of a DNA sequence. |
| `gc_content.py` | Calculates GC content as a percentage. |
| `dna_to_rna.py` | Transcribes DNA to RNA (T to U). |
| `reverse_dna.py` | Reverses a DNA sequence. |
| `reverse_complement.py` | Returns the reverse complement. |
| `rna_translation.py` | Translates RNA to protein with 1-letter, 3-letter and full amino acid names. |

Each takes an optional sequence argument and falls back to a short example:

```bash
python gc_content.py
python gc_content.py ATCGGCTAGCGCGATAT
```
