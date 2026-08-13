#!/usr/bin/env python3
"""
Translate an RNA sequence into protein.

Prints each codon with its 1-letter, 3-letter and full amino acid name, and
stops at the first stop codon.
"""

import argparse

CODONS = {
    "UUU": ("F", "Phe", "Phenylalanine"), "UUC": ("F", "Phe", "Phenylalanine"),
    "UUA": ("L", "Leu", "Leucine"), "UUG": ("L", "Leu", "Leucine"),
    "CUU": ("L", "Leu", "Leucine"), "CUC": ("L", "Leu", "Leucine"),
    "CUA": ("L", "Leu", "Leucine"), "CUG": ("L", "Leu", "Leucine"),
    "AUU": ("I", "Ile", "Isoleucine"), "AUC": ("I", "Ile", "Isoleucine"),
    "AUA": ("I", "Ile", "Isoleucine"),
    "AUG": ("M", "Met", "Methionine"),
    "GUU": ("V", "Val", "Valine"), "GUC": ("V", "Val", "Valine"),
    "GUA": ("V", "Val", "Valine"), "GUG": ("V", "Val", "Valine"),
    "UCU": ("S", "Ser", "Serine"), "UCC": ("S", "Ser", "Serine"),
    "UCA": ("S", "Ser", "Serine"), "UCG": ("S", "Ser", "Serine"),
    "AGU": ("S", "Ser", "Serine"), "AGC": ("S", "Ser", "Serine"),
    "CCU": ("P", "Pro", "Proline"), "CCC": ("P", "Pro", "Proline"),
    "CCA": ("P", "Pro", "Proline"), "CCG": ("P", "Pro", "Proline"),
    "ACU": ("T", "Thr", "Threonine"), "ACC": ("T", "Thr", "Threonine"),
    "ACA": ("T", "Thr", "Threonine"), "ACG": ("T", "Thr", "Threonine"),
    "GCU": ("A", "Ala", "Alanine"), "GCC": ("A", "Ala", "Alanine"),
    "GCA": ("A", "Ala", "Alanine"), "GCG": ("A", "Ala", "Alanine"),
    "UAU": ("Y", "Tyr", "Tyrosine"), "UAC": ("Y", "Tyr", "Tyrosine"),
    "CAU": ("H", "His", "Histidine"), "CAC": ("H", "His", "Histidine"),
    "CAA": ("Q", "Gln", "Glutamine"), "CAG": ("Q", "Gln", "Glutamine"),
    "AAU": ("N", "Asn", "Asparagine"), "AAC": ("N", "Asn", "Asparagine"),
    "AAA": ("K", "Lys", "Lysine"), "AAG": ("K", "Lys", "Lysine"),
    "GAU": ("D", "Asp", "Aspartic acid"), "GAC": ("D", "Asp", "Aspartic acid"),
    "GAA": ("E", "Glu", "Glutamic acid"), "GAG": ("E", "Glu", "Glutamic acid"),
    "UGU": ("C", "Cys", "Cysteine"), "UGC": ("C", "Cys", "Cysteine"),
    "UGG": ("W", "Trp", "Tryptophan"),
    "CGU": ("R", "Arg", "Arginine"), "CGC": ("R", "Arg", "Arginine"),
    "CGA": ("R", "Arg", "Arginine"), "CGG": ("R", "Arg", "Arginine"),
    "AGA": ("R", "Arg", "Arginine"), "AGG": ("R", "Arg", "Arginine"),
    "GGU": ("G", "Gly", "Glycine"), "GGC": ("G", "Gly", "Glycine"),
    "GGA": ("G", "Gly", "Glycine"), "GGG": ("G", "Gly", "Glycine"),
    "UAA": ("*", "Stop", "Stop"), "UAG": ("*", "Stop", "Stop"),
    "UGA": ("*", "Stop", "Stop"),
}


def translate_rna(rna, verbose=True):
    """Translate RNA to protein. Returns (1-letter, 3-letter, full-name) strings."""
    one_letter, three_letter, full_names = [], [], []

    if verbose:
        print("RNA sequence:", rna)
        print()

    for position in range(0, len(rna), 3):
        codon = rna[position:position + 3]
        if len(codon) < 3:
            break

        if codon not in CODONS:
            if verbose:
                print(codon, "-> unknown codon")
            continue

        code_1, code_3, full = CODONS[codon]
        if verbose:
            print("{} -> {} ({}) - {}".format(codon, code_1, code_3, full))

        if code_1 == "*":
            break

        one_letter.append(code_1)
        three_letter.append(code_3)
        full_names.append(full)

    return "".join(one_letter), "-".join(three_letter), " - ".join(full_names)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rna", nargs="?", default="AUGGCUUUU",
                        help="RNA sequence (default: a short example)")
    args = parser.parse_args()

    one, three, full = translate_rna(args.rna.upper())

    print()
    print("Protein (1-letter):", one)
    print("Protein (3-letter):", three)
    print("Protein (full names):", full)


if __name__ == "__main__":
    main()
