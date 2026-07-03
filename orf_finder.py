# ORF Finder - FASTA Notebook Version (updated)
# This version is designed for Jupyter Notebook or plain Python.
# It does NOT use input() or argparse, because those can get stuck in notebooks.
# It analyzes all sequences in a FASTA file and searches both forward strand
# and reverse complement.
#
# Improvements over the basic version:
#   1. Minimum ORF length filter (in amino acids), applied at the source.
#   2. Reverse-strand coordinates are mapped back to the ORIGINAL sequence,
#      so start/end always refer to the numbering in your input FASTA.
#   3. O(n) reverse complement using str.translate instead of O(n^2) string
#      concatenation.

import csv
import os


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


STOP_CODONS = ["UAA", "UAG", "UGA"]


AMINO_ACID_WEIGHTS = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10,
    "C": 121.16, "Q": 146.15, "E": 147.13, "G": 75.07,
    "H": 155.16, "I": 131.17, "L": 131.17, "K": 146.19,
    "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09,
    "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}


# Minimum number of amino acids for an ORF to be reported.
# A common rule of thumb is 100 aa (~300 bp) to filter out random ORFs.
# Set to 0 to keep every ORF (useful for the tiny example sequences below).
MIN_PROTEIN_LENGTH_AA = 100


# Translation table for the fast reverse complement (built once).
# Extend both strings if you later want to allow IUPAC codes, e.g. "ACGTN" -> "TGCAN".
_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")


def clean_dna_sequence(dna_sequence):
    """Clean DNA sequence text by removing spaces, newlines and tabs."""
    dna_sequence = dna_sequence.upper()
    dna_sequence = dna_sequence.replace(" ", "")
    dna_sequence = dna_sequence.replace("\n", "")
    dna_sequence = dna_sequence.replace("\t", "")
    return dna_sequence


def is_valid_dna(dna_sequence):
    """Return True if the sequence only contains A, T, C and G."""
    if len(dna_sequence) == 0:
        return False

    for base in dna_sequence:
        if base not in ["A", "T", "C", "G"]:
            return False

    return True


def dna_to_rna(dna_sequence):
    """Convert DNA to RNA by replacing thymine with uracil."""
    return dna_sequence.replace("T", "U")


def reverse_complement(dna_sequence):
    """
    Create the reverse complement DNA strand in O(n).

    str.translate does the base complement in a single C-level pass, and [::-1]
    reverses the result. This replaces the O(n^2) character-by-character string
    concatenation used in the original version.
    """
    return dna_sequence.translate(_COMPLEMENT_TABLE)[::-1]


def gc_content(dna_sequence):
    """Calculate GC content percentage."""
    if len(dna_sequence) == 0:
        return 0

    g = dna_sequence.count("G")
    c = dna_sequence.count("C")

    return (g + c) / len(dna_sequence) * 100


def translate_rna_orf(rna_orf):
    """Translate an RNA ORF into 1-letter and 3-letter protein codes."""
    protein_one_letter = ""
    protein_three_letter_list = []

    for position in range(0, len(rna_orf), 3):
        codon = rna_orf[position:position + 3]

        if codon not in CODON_TABLE:
            continue

        amino_acid = CODON_TABLE[codon]
        one_letter_code = amino_acid[0]
        three_letter_code = amino_acid[1]

        if one_letter_code == "*":
            break

        protein_one_letter += one_letter_code
        protein_three_letter_list.append(three_letter_code)

    protein_three_letter = "-".join(protein_three_letter_list)

    return protein_one_letter, protein_three_letter


def protein_molecular_weight(protein_sequence):
    """
    Estimate protein molecular weight in Daltons.

    The amino acid table contains approximate free amino acid masses.
    A water molecule is removed for each peptide bond during protein formation.
    """
    if len(protein_sequence) == 0:
        return 0

    weight = 0

    for amino_acid in protein_sequence:
        weight += AMINO_ACID_WEIGHTS[amino_acid]

    water_loss = 18.015 * (len(protein_sequence) - 1)

    return weight - water_loss


def find_orfs_in_dna_strand(dna_sequence, strand_name, sequence_id,
                            sequence_header, min_protein_length=0,
                            original_length=None):
    """
    Find complete ORFs in all three reading frames of one DNA strand.

    dna_sequence is the strand actually being searched. For the reverse
    complement strand, pass the reverse-complemented sequence here and give
    original_length (the length of the ORIGINAL forward sequence) so that the
    reported start/end can be mapped back to original-sequence coordinates.

    ORFs shorter than min_protein_length amino acids are skipped at the source,
    before translation and statistics, to save work on large sequences.
    """
    rna_sequence = dna_to_rna(dna_sequence)
    is_reverse = original_length is not None
    orfs = []

    for frame in range(3):
        for start_position in range(frame, len(rna_sequence) - 2, 3):
            start_codon = rna_sequence[start_position:start_position + 3]

            if start_codon == "AUG":
                for stop_position in range(start_position + 3, len(rna_sequence) - 2, 3):
                    stop_codon = rna_sequence[stop_position:stop_position + 3]

                    if stop_codon in STOP_CODONS:
                        # Amino acids = codons between start and stop (stop excluded).
                        protein_length = (stop_position - start_position) // 3

                        if protein_length < min_protein_length:
                            break  # too short - move on to the next AUG

                        rna_orf = rna_sequence[start_position:stop_position + 3]
                        dna_orf = dna_sequence[start_position:stop_position + 3]
                        protein_one, protein_three = translate_rna_orf(rna_orf)

                        # 1-based coordinates on the strand we searched.
                        searched_start = start_position + 1
                        searched_end = stop_position + 3

                        if is_reverse:
                            # Map back to original forward-sequence coordinates.
                            # A position p on the reverse complement corresponds to
                            # (L - p + 1) on the original, and the interval flips.
                            report_start = original_length - searched_end + 1
                            report_end = original_length - searched_start + 1
                        else:
                            report_start = searched_start
                            report_end = searched_end

                        orf = {
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
                        }

                        orfs.append(orf)
                        break

    return orfs


def find_orfs_both_strands(dna_sequence, sequence_id, sequence_header,
                           min_protein_length=0):
    """Find ORFs on both the forward strand and reverse complement strand."""
    forward_orfs = find_orfs_in_dna_strand(
        dna_sequence,
        "Forward",
        sequence_id,
        sequence_header,
        min_protein_length,
    )

    reverse_complement_dna = reverse_complement(dna_sequence)

    reverse_orfs = find_orfs_in_dna_strand(
        reverse_complement_dna,
        "Reverse complement",
        sequence_id,
        sequence_header,
        min_protein_length,
        original_length=len(dna_sequence),
    )

    return forward_orfs + reverse_orfs


def add_orf_statistics(orfs):
    """Add length, GC content, molecular weight and protein length fields."""
    if len(orfs) == 0:
        return orfs

    for orf in orfs:
        dna_orf = orf["dna"]
        protein = orf["protein_1_letter"]

        orf["orf_length_bp"] = len(dna_orf)
        orf["protein_length_aa"] = len(protein)
        orf["gc_content"] = round(gc_content(dna_orf), 2)
        orf["protein_molecular_weight"] = round(protein_molecular_weight(protein), 2)
        orf["longest_orf"] = False

    return orfs


def collapse_nested_orfs(orfs):
    """
    Keep only the longest ORF per stop codon.

    ORFs that share the same stop codon are the same gene read from different
    internal start codons (AUG). The longest one starts at the first AUG; the
    shorter ones are nested copies of the same protein. This keeps the longest
    and drops the rest, so each real gene is reported once - the way real ORF
    finders report one ORF per stop codon.
    """
    longest_by_stop = {}

    for orf in orfs:
        # The stop codon sits at the high-coordinate end on the forward strand.
        # For the reverse strand the coordinates were mapped back to the original
        # sequence, so the stop ends up at the low-coordinate end instead.
        if orf["strand"] == "Forward":
            stop_coordinate = orf["end"]
        else:
            stop_coordinate = orf["start"]

        key = (orf["sequence_id"], orf["strand"], stop_coordinate)
        length = orf["end"] - orf["start"] + 1

        current = longest_by_stop.get(key)
        if current is None or length > (current["end"] - current["start"] + 1):
            longest_by_stop[key] = orf

    return list(longest_by_stop.values())


def add_longest_orf_flags_per_sequence(orfs):
    """Mark the longest ORF separately for each sequence in the FASTA file."""
    if len(orfs) == 0:
        return orfs

    orfs_by_sequence = {}

    for orf in orfs:
        sequence_id = orf["sequence_id"]

        if sequence_id not in orfs_by_sequence:
            orfs_by_sequence[sequence_id] = []

        orfs_by_sequence[sequence_id].append(orf)

    for sequence_orfs in orfs_by_sequence.values():
        sequence_orfs.sort(key=lambda orf: orf["orf_length_bp"], reverse=True)

        for index, orf in enumerate(sequence_orfs):
            orf["longest_orf"] = index == 0

    return orfs


def parse_fasta(filepath):
    """Read a FASTA file and return a list of (sequence_id, header, sequence) tuples."""
    sequences = []
    current_header = None
    current_sequence = ""

    with open(filepath, "r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()

            if line == "":
                continue

            if line.startswith(">"):
                if current_header is not None:
                    sequence = clean_dna_sequence(current_sequence)
                    sequence_id = current_header.split()[0]
                    sequences.append((sequence_id, current_header, sequence))

                current_header = line[1:].strip()
                current_sequence = ""

            else:
                current_sequence += line

    if current_header is not None:
        sequence = clean_dna_sequence(current_sequence)
        sequence_id = current_header.split()[0]
        sequences.append((sequence_id, current_header, sequence))

    return sequences


def analyze_fasta_file(filepath, min_protein_length=MIN_PROTEIN_LENGTH_AA):
    """Analyze all sequences in a FASTA file and return ORF results."""
    sequences = parse_fasta(filepath)

    if len(sequences) == 0:
        print("No sequences found in FASTA file:", filepath)
        return []

    all_orfs = []
    skipped_sequences = []

    print("=" * 50)
    print("FASTA file:", filepath)
    print("Sequences found:", len(sequences))
    print("Minimum ORF length:", min_protein_length, "aa")
    print("=" * 50)

    for sequence_id, sequence_header, dna_sequence in sequences:
        print()
        print("Sequence:", sequence_id)
        print("Header:", sequence_header)
        print("Length:", len(dna_sequence), "bp")

        if not is_valid_dna(dna_sequence):
            print("Skipped: invalid DNA sequence. Please use only A, T, C, and G.")
            skipped_sequences.append(sequence_id)
            continue

        orfs = find_orfs_both_strands(
            dna_sequence, sequence_id, sequence_header, min_protein_length
        )
        raw_count = len(orfs)
        orfs = collapse_nested_orfs(orfs)
        orfs = add_orf_statistics(orfs)
        all_orfs.extend(orfs)

        print("ORFs found in this sequence:", len(orfs),
              "(collapsed from", raw_count, "nested ORFs)")

    if len(skipped_sequences) > 0:
        print()
        print("Skipped invalid sequences:", ", ".join(skipped_sequences))

    if len(all_orfs) > 0:
        all_orfs = add_longest_orf_flags_per_sequence(all_orfs)

    return all_orfs


def summarize_orfs(orfs):
    """Print summary of all ORF results."""
    print()
    print("=" * 50)
    print("ORF Finder Results")
    print("=" * 50)

    forward_count = 0
    reverse_count = 0
    sequence_ids = set()

    for orf in orfs:
        sequence_ids.add(orf["sequence_id"])

        if orf["strand"] == "Forward":
            forward_count += 1
        else:
            reverse_count += 1

    print("Sequences analyzed:", len(sequence_ids))
    print("Total ORFs found:", len(orfs))
    print("Forward strand ORFs:", forward_count)
    print("Reverse complement ORFs:", reverse_count)
    print("=" * 50)


def print_orfs(orfs):
    """Print ORF results in a readable format."""
    if len(orfs) == 0:
        print("No complete ORFs found (try lowering the minimum ORF length).")
        return

    summarize_orfs(orfs)

    for orf_number, orf in enumerate(orfs, start=1):
        print()
        print("ORF", orf_number)

        if orf["longest_orf"]:
            print("Longest ORF in sequence", orf["sequence_id"])

        print("-" * 40)
        print("Sequence ID:", orf["sequence_id"])
        print("Header:", orf["sequence_header"])
        print("Strand:", orf["strand"])
        print("Reading frame:", orf["frame"])
        # start/end are always in ORIGINAL sequence coordinates now.
        print("Position on original sequence:", str(orf["start"]) + "-" + str(orf["end"]))
        print("DNA ORF (coding strand 5'->3'):", orf["dna"])
        print("RNA ORF:", orf["rna"])
        print("ORF length:", orf["orf_length_bp"], "bp")
        print("GC content:", orf["gc_content"], "%")
        print("Protein 1-letter:", orf["protein_1_letter"])
        print("Protein 3-letter:", orf["protein_3_letter"])
        print("Protein length:", orf["protein_length_aa"], "aa")
        print("Protein molecular weight:", orf["protein_molecular_weight"], "Da")


def export_orfs_to_csv(orfs, filename="orf_results.csv"):
    """Export ORF results to a CSV file."""
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "sequence_id",
            "sequence_header",
            "strand",
            "frame",
            "start",
            "end",
            "orf_length_bp",
            "gc_content",
            "protein_1_letter",
            "protein_3_letter",
            "protein_length_aa",
            "protein_molecular_weight",
            "longest_orf",
            "dna",
            "rna",
        ]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for orf in orfs:
            writer.writerow(orf)

    print("Results exported to", filename)


def build_output_csv_path(fasta_filepath, output_csv=None):
    """Create a default CSV filename if the user has not provided one."""
    if output_csv is not None:
        return output_csv

    base_name = os.path.splitext(os.path.basename(fasta_filepath))[0]

    return base_name + "_orf_results.csv"


def create_example_fasta(filename="example_sequences.fasta"):
    """Create a small example FASTA file so the notebook can be tested immediately."""
    with open(filename, "w", encoding="utf-8") as fasta_file:
        fasta_file.write(">example_sequence_1 multiple ORFs\n")
        fasta_file.write("GGGCCCATGGCTTTTGGAGAATAACCCGGGATGAAACCCGGGTAGTTTCCCATGTTTGGATGAGGGAAATGCCCGGCTTCTAA\n")
        fasta_file.write(">example_sequence_2 short test sequence\n")
        fasta_file.write("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n")

    print("Example FASTA file created:", filename)


def run_orf_finder(fasta_filepath, output_csv=None,
                   min_protein_length=MIN_PROTEIN_LENGTH_AA):
    """Run the FASTA ORF Finder from a notebook-friendly function."""
    if fasta_filepath is None or fasta_filepath == "":
        print("No FASTA file provided.")
        return []

    if not os.path.isfile(fasta_filepath):
        print("File not found:", fasta_filepath)
        return []

    orfs = analyze_fasta_file(fasta_filepath, min_protein_length)
    print_orfs(orfs)

    if len(orfs) > 0:
        output_csv_path = build_output_csv_path(fasta_filepath, output_csv)
        export_orfs_to_csv(orfs, output_csv_path)

    return orfs


# ------------------------------------------------------------
# Notebook run section
# ------------------------------------------------------------
# Leave USE_EXAMPLE_FASTA = True if you want to test the program immediately.
# When you want to analyze your own FASTA file:
# 1. Put the FASTA file in the same folder as this notebook
# 2. Change USE_EXAMPLE_FASTA to False
# 3. Change FASTA_FILE to your filename, for example "sars_cov_2.fasta"
#
# Note: the example sequences are tiny, so with MIN_PROTEIN_LENGTH_AA = 100 they
# produce zero ORFs. For the demo we lower the threshold via MIN_ORF_LENGTH below.
# For real genomes, use MIN_PROTEIN_LENGTH_AA (100) or your own value.

# This block only runs when you execute this file DIRECTLY
# (python orf_finder.py). It does NOT run when the file is imported from a
# notebook, so importing no longer triggers the example demo.
if __name__ == "__main__":
    USE_EXAMPLE_FASTA = True
    FASTA_FILE = "example_sequences.fasta"
    OUTPUT_CSV = None
    MIN_ORF_LENGTH = 1  # aa; raise to 100 for real data

    if USE_EXAMPLE_FASTA:
        create_example_fasta(FASTA_FILE)

    orfs = run_orf_finder(FASTA_FILE, OUTPUT_CSV, min_protein_length=MIN_ORF_LENGTH)
