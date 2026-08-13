#!/usr/bin/env python3
"""Count how many times each base appears in a DNA sequence."""

import argparse


def base_counter(dna):
    counts = {base: dna.count(base) for base in "ATCG"}
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dna", nargs="?", default="ATCGGCTAGCGCGATAT",
                        help="DNA sequence (default: a short example)")
    args = parser.parse_args()
    dna = args.dna.upper()

    print("DNA sequence:", dna)
    print("Total bases:", len(dna))
    print()
    for base, count in base_counter(dna).items():
        print(base, "=", count)


if __name__ == "__main__":
    main()
