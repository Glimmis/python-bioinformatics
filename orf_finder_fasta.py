# ORF Finder - FASTA Notebook Version
# This version is designed for Jupyter Notebook.
# It does NOT use input() or argparse, because those can get stuck in notebooks.
# It analyzes all sequences in a FASTA file and searches both forward strand and reverse complement.

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
    """Create the reverse complement DNA strand."""
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }

    reversed_complement = ""

    for base in dna_sequence:
        reversed_complement = complement[base] + reversed_complement

    return reversed_complement


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


def find_orfs_in_dna_strand(dna_sequence, strand_name, sequence_id, sequence_header):
    """Find complete ORFs in all three reading frames of one DNA strand."""
    rna_sequence = dna_to_rna(dna_sequence)
    orfs = []

    for frame in range(3):
        for start_position in range(frame, len(rna_sequence) - 2, 3):
            start_codon = rna_sequence[start_position:start_position + 3]

            if start_codon == "AUG":
                for stop_position in range(start_position + 3, len(rna_sequence) - 2, 3):
                    stop_codon = rna_sequence[stop_position:stop_position + 3]

                    if stop_codon in STOP_CODONS:
                        rna_orf = rna_sequence[start_position:stop_position + 3]
                        dna_orf = dna_sequence[start_position:stop_position + 3]
                        protein_one, protein_three = translate_rna_orf(rna_orf)

                        orf = {
                            "sequence_id": sequence_id,
                            "sequence_header": sequence_header,
                            "strand": strand_name,
                            "frame": frame + 1,
                            "start": start_position + 1,
                            "end": stop_position + 3,
                            "dna": dna_orf,
                            "rna": rna_orf,
                            "protein_1_letter": protein_one,
                            "protein_3_letter": protein_three,
                        }

                        orfs.append(orf)
                        break

    return orfs


def find_orfs_both_strands(dna_sequence, sequence_id, sequence_header):
    """Find ORFs on both the forward strand and reverse complement strand."""
    forward_orfs = find_orfs_in_dna_strand(
        dna_sequence,
        "Forward",
        sequence_id,
        sequence_header
    )

    reverse_complement_dna = reverse_complement(dna_sequence)

    reverse_orfs = find_orfs_in_dna_strand(
        reverse_complement_dna,
        "Reverse complement",
        sequence_id,
        sequence_header
    )

    return forward_orfs + reverse_orfs


def add_orf_statistics(orfs):
    """Add length, GC content, molecular weight and longest-ORF flag."""
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

    orfs.sort(key=lambda orf: orf["orf_length_bp"], reverse=True)
    orfs[0]["longest_orf"] = True

    return orfs


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


def analyze_fasta_file(filepath):
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

        orfs = find_orfs_both_strands(dna_sequence, sequence_id, sequence_header)
        orfs = add_orf_statistics(orfs)
        all_orfs.extend(orfs)

        print("ORFs found in this sequence:", len(orfs))

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
        print("No complete ORFs found.")
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
        print("Position on searched strand:", str(orf["start"]) + "-" + str(orf["end"]))
        print("DNA ORF:", orf["dna"])
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


def run_orf_finder(fasta_filepath, output_csv=None):
    """Run the FASTA ORF Finder from a notebook-friendly function."""
    if fasta_filepath is None or fasta_filepath == "":
        print("No FASTA file provided.")
        return []

    if not os.path.isfile(fasta_filepath):
        print("File not found:", fasta_filepath)
        return []

    orfs = analyze_fasta_file(fasta_filepath)
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

USE_EXAMPLE_FASTA = True
FASTA_FILE = "example_sequences.fasta"
OUTPUT_CSV = None

if USE_EXAMPLE_FASTA:
    create_example_fasta(FASTA_FILE)

orfs = run_orf_finder(FASTA_FILE, OUTPUT_CSV)
