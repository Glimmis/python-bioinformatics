#!/usr/bin/env python3
"""Return the reverse complement of a DNA sequence."""

import argparse

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(dna):
    """O(n) via str.translate, rather than character-by-character concatenation."""
    return dna.translate(COMPLEMENT)[::-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dna", nargs="?", default="ATCGGCTAGCGCGATAT",
                        help="DNA sequence (default: a short example)")
    args = parser.parse_args()
    dna = args.dna.upper()

    print("Original:          ", dna)
    print("Reverse complement:", reverse_complement(dna))


if __name__ == "__main__":
    main()
