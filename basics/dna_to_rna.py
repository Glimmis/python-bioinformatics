#!/usr/bin/env python3
"""Transcribe DNA to RNA by replacing every T with U."""

import argparse


def dna_to_rna(dna):
    return dna.replace("T", "U")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dna", nargs="?", default="ATCGGCTAGCGCGATAT",
                        help="DNA sequence (default: a short example)")
    args = parser.parse_args()
    dna = args.dna.upper()

    print("DNA:", dna)
    print("RNA:", dna_to_rna(dna))


if __name__ == "__main__":
    main()
