#!/usr/bin/env python3
"""Return the length of a DNA sequence in bases."""

import argparse


def dna_length(dna):
    return len(dna)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dna", nargs="?", default="ATCGGCTAGCGCGATAT",
                        help="DNA sequence (default: a short example)")
    args = parser.parse_args()
    dna = args.dna.upper()

    print("DNA sequence:", dna)
    print("Length:", dna_length(dna))


if __name__ == "__main__":
    main()
