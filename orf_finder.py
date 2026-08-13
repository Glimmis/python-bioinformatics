#!/usr/bin/env python3
"""
ORF finder: structural annotation of a FASTA file.

Scans every sequence in a FASTA file for complete open reading frames (ATG to
stop codon) in all three reading frames on both the forward strand and the
reverse complement, translates them, and writes the results to CSV.

Design notes:
  * Reverse-strand coordinates are mapped back to the ORIGINAL sequence, so
    start/end always refer to the numbering in the input FASTA.
  * ORFs sharing a stop codon are collapsed to the longest one, so each gene is
    reported once - the way real ORF finders behave.
  * The minimum length filter is applied before translation, so large genomes
    do not pay for ORFs that will be discarded.
  * Reverse complement uses str.translate (O(n)) rather than character-by-
    character concatenation (O(n^2)).

Usage:
    python orf_finder.py examples/example_contig.fasta
    python orf_finder.py genome.fasta --min-length 100 -o genome_orfs.csv

Also importable:
    from orf_finder import run_orf_finder
    orfs = run_orf_finder("genome.fasta", min_protein_length=100)
"""

import argparse
import csv
import os
import sys

CODON_TABLE = {
    "UUU": ["F", "Phe"], "UUC": ["F", "Phe"],
    "UUA": ["L", "Leu"], "UUG": ["L", "Leu"],
    "UCU": ["S", "Ser"], "UCC": ["S", "Ser"], "UCA": ["S", "Ser"], "UCG": ["S", "Ser"],
    "UAU": ["Y", "Tyr"], "UAC": ["Y", "Tyr"],
    "UAA": ["*", "Stop"], "UAG": ["*", "Stop"],
    "UGU": ["C", "Cys"], "UGC": ["C", "Cys"],
    "UGA": ["*", "Stop"], "UGG": ["W", "Trp"],
    "CUU": ["L", "Leu"], "CUC": ["L", "Leu"], "CUA": ["L", "Leu"], "CUG": ["L", "Leu"],
    "CCU": ["P", "Pro"], "CCC": ["P", "Pro"], "CCA": ["P", "Pro"], "CCG": ["P", "Pro"],
    "CAU": ["H", "His"], "CAC": ["H", "His"],
    "CAA": ["Q", "Gln"], "CAG": ["Q", "Gln"],
    "CGU": ["R", "Arg"], "CGC": ["R", "Arg"], "CGA": ["R", "Arg"], "CGG": ["R", "Arg"],
    "AUU": ["I", "Ile"], "AUC": ["I", "Ile"], "AUA": ["I", "Ile"],
    "AUG": ["M", "Met"],
    "ACU": ["T", "Thr"], "ACC": ["T", "Thr"], "ACA": ["T", "Thr"], "ACG": ["T", "Thr"],
    "AAU": ["N", "Asn"], "AAC": ["N", "Asn"],
    "AAA": ["K", "Lys"], "AAG": ["K", "Lys"],
    "AGU": ["S", "Ser"], "AGC": ["S", "Ser"],
    "AGA": ["R", "Arg"], "AGG": ["R", "Arg"],
    "GUU": ["V", "Val"], "GUC": ["V", "Val"], "GUA": ["V", "Val"], "GUG": ["V", "Val"],
    "GCU": ["A", "Ala"], "GCC": ["A", "Ala"], "GCA": ["A", "Ala"], "GCG": ["A", "Ala"],
    "GAU": ["D", "Asp"], "GAC": ["D", "Asp"],
    "GAA": ["E", "Glu"], "GAG": ["E", "Glu"],
    "GGU": ["G", "Gly"], "GGC": ["G", "Gly"], "GGA": ["G", "Gly"], "GGG": ["G", "Gly"],
}

STOP_CODONS = ("UAA", "UAG", "UGA")

AMINO_ACID_WEIGHTS = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10,
    "C": 121.16, "Q": 146.15, "E": 147.13, "G": 75.07,
    "H": 155.16, "I": 131.17, "L": 131.17, "K": 146.19,
    "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09,
    "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}

# Default: 100 aa (~300 bp) is the usual rule of thumb for filtering out ORFs
# that occur by chance in random sequence.
DEFAULT_MIN_PROTEIN_LENGTH_AA = 100

_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")

CSV_FIELDS = [
    "sequence_id", "sequence_header", "strand", "frame", "start", "end",
    "orf_length_bp", "gc_content", "protein_1_letter", "protein_3_letter",
    "protein_length_aa", "protein_molecular_weight", "longest_orf", "dna", "rna",
]


def clean_dna_sequence(dna_sequence):
    """Strip whitespace and upper-case a DNA sequence."""
    return dna_sequence.upper().replace(" ", "").replace("\n", "").replace("\t", "")


def is_valid_dna(dna_sequence):
    """Return True if the sequence is non-empty and contains only A, T, C, G."""
    return len(dna_sequence) > 0 and set(dna_sequence) <= {"A", "T", "C", "G"}


def dna_to_rna(dna_sequence):
    """Transcribe DNA to RNA (T -> U)."""
    return dna_sequence.replace("T", "U")


def reverse_complement(dna_sequence):
    """Reverse complement in O(n) using str.translate."""
    return dna_sequence.translate(_COMPLEMENT_TABLE)[::-1]


def gc_content(dna_sequence):
    """GC content as a percentage."""
    if not dna_sequence:
        return 0.0
    return (dna_sequence.count("G") + dna_sequence.count("C")) / len(dna_sequence) * 100


def translate_rna_orf(rna_orf):
    """Translate an RNA ORF, returning (1-letter, 3-letter) protein strings."""
    one_letter = ""
    three_letter = []

    for position in range(0, len(rna_orf), 3):
        codon = rna_orf[position:position + 3]
        if codon not in CODON_TABLE:
            continue
        code_1, code_3 = CODON_TABLE[codon]
        if code_1 == "*":
            break
        one_letter += code_1
        three_letter.append(code_3)

    return one_letter, "-".join(three_letter)


def protein_molecular_weight(protein_sequence):
    """Estimate molecular weight in Da (free amino acid masses minus peptide water)."""
    if not protein_sequence:
        return 0.0
    weight = sum(AMINO_ACID_WEIGHTS[aa] for aa in protein_sequence)
    return weight - 18.015 * (len(protein_sequence) - 1)


def find_orfs_in_dna_strand(dna_sequence, strand_name, sequence_id, sequence_header,
                            min_protein_length=0, original_length=None):
    """
    Find complete ORFs in all three reading frames of one strand.

    Pass the reverse-complemented sequence plus original_length (length of the
    forward sequence) to have coordinates mapped back to forward numbering.
    """
    rna_sequence = dna_to_rna(dna_sequence)
    is_reverse = original_length is not None
    orfs = []

    for frame in range(3):
        for start_position in range(frame, len(rna_sequence) - 2, 3):
            if rna_sequence[start_position:start_position + 3] != "AUG":
                continue

            for stop_position in range(start_position + 3, len(rna_sequence) - 2, 3):
                if rna_sequence[stop_position:stop_position + 3] not in STOP_CODONS:
                    continue

                protein_length = (stop_position - start_position) // 3
                if protein_length < min_protein_length:
                    break

                rna_orf = rna_sequence[start_position:stop_position + 3]
                dna_orf = dna_sequence[start_position:stop_position + 3]
                protein_one, protein_three = translate_rna_orf(rna_orf)

                searched_start = start_position + 1
                searched_end = stop_position + 3

                if is_reverse:
                    report_start = original_length - searched_end + 1
                    report_end = original_length - searched_start + 1
                else:
                    report_start = searched_start
                    report_end = searched_end

                orfs.append({
                    "sequence_id": sequence_id,
                    "sequence_header": sequence_header,
                    "strand": strand_name,
                    "frame": frame + 1,
                    "start": report_start,
                    "end": report_end,
                    "dna": dna_orf,
                    "rna": rna_orf,
                    "protein_1_letter": protein_one,
                    "protein_3_letter": protein_three,
                })
                break

    return orfs


def find_orfs_both_strands(dna_sequence, sequence_id, sequence_header,
                           min_protein_length=0):
    """Search forward strand and reverse complement."""
    forward = find_orfs_in_dna_strand(
        dna_sequence, "Forward", sequence_id, sequence_header, min_protein_length)
    reverse = find_orfs_in_dna_strand(
        reverse_complement(dna_sequence), "Reverse complement", sequence_id,
        sequence_header, min_protein_length, original_length=len(dna_sequence))
    return forward + reverse


def add_orf_statistics(orfs):
    """Add length, GC content, protein length and molecular weight."""
    for orf in orfs:
        protein = orf["protein_1_letter"]
        orf["orf_length_bp"] = len(orf["dna"])
        orf["protein_length_aa"] = len(protein)
        orf["gc_content"] = round(gc_content(orf["dna"]), 2)
        orf["protein_molecular_weight"] = round(protein_molecular_weight(protein), 2)
        orf["longest_orf"] = False
    return orfs


def collapse_nested_orfs(orfs):
    """
    Keep only the longest ORF per stop codon.

    ORFs sharing a stop codon are the same gene read from different internal
    start codons. The longest starts at the first ATG; the rest are nested
    copies of the same protein.
    """
    longest_by_stop = {}

    for orf in orfs:
        # On the forward strand the stop sits at the high coordinate; after
        # mapping reverse coordinates back, it sits at the low coordinate.
        stop_coordinate = orf["end"] if orf["strand"] == "Forward" else orf["start"]
        key = (orf["sequence_id"], orf["strand"], stop_coordinate)
        length = orf["end"] - orf["start"] + 1

        current = longest_by_stop.get(key)
        if current is None or length > (current["end"] - current["start"] + 1):
            longest_by_stop[key] = orf

    return list(longest_by_stop.values())


def add_longest_orf_flags_per_sequence(orfs):
    """Flag the longest ORF within each input sequence."""
    by_sequence = {}
    for orf in orfs:
        by_sequence.setdefault(orf["sequence_id"], []).append(orf)

    for sequence_orfs in by_sequence.values():
        sequence_orfs.sort(key=lambda o: o["orf_length_bp"], reverse=True)
        for index, orf in enumerate(sequence_orfs):
            orf["longest_orf"] = index == 0

    return orfs


def parse_fasta(filepath):
    """Read a FASTA file into a list of (sequence_id, header, sequence) tuples."""
    sequences = []
    header = None
    chunks = []

    def flush():
        if header is not None:
            sequence = clean_dna_sequence("".join(chunks))
            sequences.append((header.split()[0], header, sequence))

    with open(filepath, "r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    flush()
    return sequences


def analyze_fasta_file(filepath, min_protein_length=DEFAULT_MIN_PROTEIN_LENGTH_AA,
                       quiet=False):
    """Analyze every sequence in a FASTA file and return the ORF list."""
    sequences = parse_fasta(filepath)

    if not sequences:
        print("No sequences found in FASTA file:", filepath, file=sys.stderr)
        return []

    all_orfs = []
    skipped = []

    if not quiet:
        print("=" * 50)
        print("FASTA file:", filepath)
        print("Sequences found:", len(sequences))
        print("Minimum ORF length:", min_protein_length, "aa")
        print("=" * 50)

    for sequence_id, header, dna_sequence in sequences:
        if not quiet:
            print()
            print("Sequence:", sequence_id)
            print("Header:", header)
            print("Length:", len(dna_sequence), "bp")

        if not is_valid_dna(dna_sequence):
            print("Skipped: invalid DNA sequence (only A, T, C, G allowed).",
                  file=sys.stderr)
            skipped.append(sequence_id)
            continue

        orfs = find_orfs_both_strands(dna_sequence, sequence_id, header,
                                      min_protein_length)
        raw_count = len(orfs)
        orfs = add_orf_statistics(collapse_nested_orfs(orfs))
        all_orfs.extend(orfs)

        if not quiet:
            print("ORFs found:", len(orfs),
                  "(collapsed from", raw_count, "nested ORFs)")

    if skipped:
        print("Skipped invalid sequences:", ", ".join(skipped), file=sys.stderr)

    if all_orfs:
        all_orfs = add_longest_orf_flags_per_sequence(all_orfs)

    return all_orfs


def summarize_orfs(orfs):
    """Print a short summary of the ORF results."""
    forward = sum(1 for o in orfs if o["strand"] == "Forward")
    sequence_ids = {o["sequence_id"] for o in orfs}

    print()
    print("=" * 50)
    print("ORF finder results")
    print("=" * 50)
    print("Sequences analyzed:", len(sequence_ids))
    print("Total ORFs found:", len(orfs))
    print("Forward strand ORFs:", forward)
    print("Reverse complement ORFs:", len(orfs) - forward)
    print("=" * 50)


def print_orfs(orfs, detail=False):
    """Print results; with detail=True print every ORF in full."""
    if not orfs:
        print("No complete ORFs found (try lowering --min-length).")
        return

    summarize_orfs(orfs)

    if not detail:
        return

    for number, orf in enumerate(orfs, start=1):
        print()
        print("ORF", number)
        if orf["longest_orf"]:
            print("Longest ORF in sequence", orf["sequence_id"])
        print("-" * 40)
        print("Sequence ID:", orf["sequence_id"])
        print("Strand:", orf["strand"])
        print("Reading frame:", orf["frame"])
        print("Position on original sequence:",
              "{}-{}".format(orf["start"], orf["end"]))
        print("ORF length:", orf["orf_length_bp"], "bp")
        print("GC content:", orf["gc_content"], "%")
        print("Protein (1-letter):", orf["protein_1_letter"])
        print("Protein length:", orf["protein_length_aa"], "aa")
        print("Protein molecular weight:", orf["protein_molecular_weight"], "Da")


def export_orfs_to_csv(orfs, filename="orf_results.csv"):
    """Write the ORF table to CSV."""
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(orfs)
    print("Results written to", filename)


def build_output_csv_path(fasta_filepath, output_csv=None):
    """Derive a default CSV filename from the FASTA filename."""
    if output_csv:
        return output_csv
    base = os.path.splitext(os.path.basename(fasta_filepath))[0]
    return base + "_orfs.csv"


def run_orf_finder(fasta_filepath, output_csv=None,
                   min_protein_length=DEFAULT_MIN_PROTEIN_LENGTH_AA,
                   write_csv=True, detail=False, quiet=False):
    """Run the full ORF search. Returns the list of ORF dicts."""
    if not fasta_filepath or not os.path.isfile(fasta_filepath):
        print("File not found:", fasta_filepath, file=sys.stderr)
        return []

    orfs = analyze_fasta_file(fasta_filepath, min_protein_length, quiet=quiet)

    if not quiet:
        print_orfs(orfs, detail=detail)

    if orfs and write_csv:
        export_orfs_to_csv(orfs, build_output_csv_path(fasta_filepath, output_csv))

    return orfs


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Find open reading frames in a FASTA file and export them to CSV.",
        epilog="Example: python orf_finder.py examples/example_contig.fasta "
               "--min-length 100 -o contig_orfs.csv",
    )
    parser.add_argument("fasta", help="input FASTA file")
    parser.add_argument("-o", "--output", default=None,
                        help="output CSV path (default: <fasta basename>_orfs.csv)")
    parser.add_argument("-m", "--min-length", type=int,
                        default=DEFAULT_MIN_PROTEIN_LENGTH_AA,
                        help="minimum protein length in amino acids (default: %(default)s)")
    parser.add_argument("--detail", action="store_true",
                        help="print every ORF in full, not just the summary")
    parser.add_argument("--no-csv", action="store_true",
                        help="print results without writing a CSV file")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="suppress progress output")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    orfs = run_orf_finder(
        args.fasta,
        output_csv=args.output,
        min_protein_length=args.min_length,
        write_csv=not args.no_csv,
        detail=args.detail,
        quiet=args.quiet,
    )
    return 0 if orfs else 1


if __name__ == "__main__":
    raise SystemExit(main())
