#!/usr/bin/env python3
"""Reverse a DNA sequence (not the complement - see reverse_complement.py)."""

import argparse


def reverse_dna(dna):
    return dna[::-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dna", nargs="?", default="ATCGGCTAGCGCGATAT",
                        help="DNA sequence (default: a short example)")
    args = parser.parse_args()
    dna = args.dna.upper()

    print("Original:", dna)
    print("Reversed:", reverse_dna(dna))


if __name__ == "__main__":
    main()
